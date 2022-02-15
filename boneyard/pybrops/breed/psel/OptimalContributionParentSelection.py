import cvxpy
import math
import numpy

from . import ParentSelectionOperator

import pybrops.core.random
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_ndarray
from pybrops.core.error import cond_check_is_ndarray
from pybrops.core.error import cond_check_is_Generator

class OptimalContributionParentSelection(ParentSelectionOperator):
    """docstring for OptimalContributionParentSelection."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, k_p, traitwt_p, inbfn_p, ncross, nprogeny, cmatcls, bvtype = "gebv", rng = None, **kwargs):
        """
        cmatcls : CoancestryMatrix class
        inbfn_p : function
            inbfn_p(t_cur, t_max)
            Returns constraint for mean population inbreeding defined as:
                (1/2) x'Ax = x'Kx <= inbfn_p(t_cur, t_max)
        bvtype : str
            Whether to use GEBVs or phenotypic EBVs.
            Options:
                "gebv"      Use GEBVs
                "ebv"       Use EBVs
        """
        super(OptimalContributionParentSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(k_p, "k_p")
        cond_check_is_ndarray(traitwt_p, "traitwt_p")
        # TODO: check inbfn_p
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        # TODO: check cmatcls
        # TODO: check bvtype
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.k_p = k_p
        self.traitwt_p = traitwt_p
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.cmatcls = cmatcls
        self.bvtype = bvtype
        self.inbfn_p = inbfn_p
        self.rng = pybrops.core.random if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def get_bv(self, geno, bval, gmod, **kwargs):
        option = self.bvtype.lower()
        if option == "gebv":
            return gmod["cand"].predict(geno["cand"]).mat   # (n,t)
        elif option == "ebv":
            return bval["cand"].mat                         # (n,t)
        else:
            raise RuntimeError("unknown bvtype")

    def calc_K(self, gmat, **kwargs):
        kmat = self.cmatcls.from_gmat(gmat)  # get coancestry (kinship) matrix
        return kmat

    def sus(self, k, contrib):
        """
        Stochastic universal sampling

        contrib : numpy.ndarray
            Contribution matrix of shape (n,).
            Where:
                'n' is the number of individuals.
            Restrictions:
                Values are restricted to [0, inf].
                Sum of values in vector must be > 0.0
        """
        tot_fit = contrib.sum()                         # calculate the total fitness
        ptr_dist = tot_fit / k                          # calculate the distance between pointers
        indices = contrib.argsort()[::-1]               # get indices for sorted individuals
        cumsum = contrib[indices].cumsum()              # get cumulative sum of elements
        offset = self.rng.uniform(0.0, ptr_dist)        # get random start point
        sel = []                                        # declare output list
        ix = 0                                          # index for cumsum
        ptrs = numpy.arange(offset, tot_fit, ptr_dist)  # create pointers
        for ptr in ptrs:                                # for each pointer
            while cumsum[ix] < ptr:                     # advance index to correct location
                ix += 1
            sel.append(indices[ix])                     # append index with element
        sel = numpy.array(sel)                          # convert to ndarray
        return sel

    def pselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        geno : dict
            A dict containing genotypic data for all breeding populations.
            Must have the following fields:
                Field | Type                         | Description
                ------+------------------------------+--------------------------
                cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                main  | PhasedGenotypeMatrix         | Main breeding population
                queue | List of PhasedGenotypeMatrix | Breeding populations on queue
                ""
        bval : dict
            A dict containing breeding value data.
            Must have the following fields:
                Field      | Type                        | Description
                -----------+-----------------------------+----------------------
                cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                main       | BreedingValueMatrix         | Main breeding population breeding values
                main_true  | BreedingValueMatrix         | Main breeding population true breeding values
        gmod : dict
            A dict containing genomic models.
            Must have the following fields:
                Field | Type                 | Description
                ------+----------------------+----------------------------------
                cand  | GenomicModel         | Parental candidate breeding population genomic model
                main  | GenomicModel         | Main breeding population genomic model
                true  | GenomicModel         | True genomic model for trait(s)
        k : int
        traitwt : numpy.ndarray
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgvmat, sel, ncross, nprogeny, misc)
            pgvmat : PhasedGenotypeVariantMatrix
                A PhasedGenotypeVariantMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgvmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_p
        if traitwt is None:
            traitwt = self.traitwt_p

        # TODO: implement multiobjective/non-dominated OCS
        if traitwt is None:
            raise RuntimeError("multi-objective optimal contribution selection not implemented")

        ##############################
        # Optimization problem setup #
        ##############################

        # get breeding values
        bv = self.get_bv(geno, bval, gmod)              # (n,t)

        # calculate kinship matrix (ndarray)
        K = self.calc_K(geno["cand"]).mat               # (n,n)

        # apply weights to get weighted BVs
        wbv = bv.dot(traitwt)                           # (n,t).(t,) -> (n,)

        # get number of taxa
        ntaxa = len(wbv)                                # scalar

        # cholesky decompose K into K = C'C
        # needed to decompose x'Kx <= inbfn_p to:
        # ||Cx||_2 <= sqrt(inbfn_p)
        # FIXME: will fail if matrix is not positive definite
        try:
            C = numpy.linalg.cholesky(K).T              # (n,n)
        except numpy.linalg.LinAlgError:
            # TODO: revert to WFS or uniform here
            raise RuntimeError(
                "Kinship matrix is not positive definite."+
                "This could be due to loss of genetic diversity."
            )

        # get sqrt(max inbreeding)
        inbmax = math.sqrt(self.inbfn_p(t_cur, t_max))  # sqrt(max inbreeding)

        # get adding vector
        ones = numpy.ones(ntaxa, dtype="float64")       # (n,)

        #######################################
        # Define and solve the CVXPY problem. #
        #######################################

        # define variable to optimize
        x = cvxpy.Variable(ntaxa)                       # (n,)

        # define the objective function
        soc_objfn = cvxpy.Maximize(wbv @ x)             # max (bv)'(sel)

        # define constraints
        soc_constraints = [
            cvxpy.SOC(inbmax, C @ x),                   # ||C @ x||_2 <= inbmax
            ones @ x == 1.0,                            # sum(x_i) == 1
            x >= 0.0                                    # x_i >= 0 for all i
        ]

        # define problem
        prob = cvxpy.Problem(
            soc_objfn,
            soc_constraints
        )

        # solve the problem
        sol = prob.solve()

        ##################
        # select parents #
        ##################

        # if the problem is unsolvable, revert to windowed fitness proportional selection
        sel = None
        if prob.status != "optimal":
            contrib = wbv - wbv.min()                   # window fitnesses
            if contrib.sum() < 1e-10:                   # if everything is near identical
                contrib = numpy.repeat(1/ntaxa, ntaxa)  # default to equal chance
            sel = self.sus(k, contrib)                  # sample indices
        else:
            sel = self.sus(k, x.value)                  # sample indices

        ################################################################
        # Make sure parents are forced to outbreed as much as possible #
        ################################################################
        self.rng.shuffle(sel)                           # start with random shuffle
        for st in range(0, len(sel), 2):                # for each female
            sp = st+1                                   # get stop index (male)
            j = st+2                                    # get one beyond stop index
            while (sel[st]==sel[sp]) and (j<len(sel)):  # if female == male && within array bounds
                sel[sp], sel[j] = sel[j], sel[sp]       # exchange with next entry
                j += 1                                  # increment beyond index

        # nothing to do... ...yet
        misc = {}

        return geno["cand"], sel, self.ncross, self.nprogeny, misc

    def pobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        bv = self.get_bv(geno, bval, gmod)              # (n,t)

        # TODO: put this in static library
        def objfn(sel, bv = bv, traitwt = traitwt, **kwargs):
            """
            Calculate genotypic mean from contributions:
                                    (bv)'(sel)

            Parameters
            ----------
            sel : numpy.ndarray
                A selection weights matrix of shape (n,)
                Where:
                    'n' is the number of individuals.
            bv : numpy.ndarray
                A breeding value matrix of shape (n, t).
                Where:
                    'n' is the number of individuals.
                    't' is the number of trait objectives.
            traitwt : numpy.ndarray, None
                A trait objective coefficients matrix of shape (t,).
                Where:
                    't' is the number of trait objectives.
                These are used to weigh objectives in the weight sum method.
                If None, do not multiply BVs by a weight sum vector.

            Returns
            -------
            avgbv : float, numpy.ndarray
                Genotypic mean for the selection.
                If traitwt is provided, returns a scalar.
                If traitwt is not provided, returns a matrix of shape (t,)
                Where:
                    't' is the number of trait objectives.
            """
            # calculate weighted average BVs
            # if traitwt is None:
            #     (n,).(n,t) -> (t,)
            # else:
            #     (n,).{(n,t).(t,) -> (n,)} -> scalar
            avgbv = sel.dot(bv) if traitwt is None else sel.dot(bv.dot(traitwt))
            return avgbv

        return objfn

    def pobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
        """
        Return a vectorized objective function.
        """
        bv = self.get_bv(geno, bval, gmod)              # (n,t)

        def objfn_vec(sel, bv = bv, traitwt = traitwt, **kwargs):
            """
            Calculate genotypic mean from contributions:
                                    (bv)'(sel)

            Parameters
            ----------
            sel : numpy.ndarray
                A selection weights matrix of shape (j,n)
                Where:
                    'j' is the number of selection configurations.
                    'n' is the number of individuals.
            bv : numpy.ndarray
                A breeding value matrix of shape (n, t).
                Where:
                    'n' is the number of individuals.
                    't' is the number of trait objectives.
            traitwt : numpy.ndarray, None
                A trait objective coefficients matrix of shape (t,).
                Where:
                    't' is the number of trait objectives.
                These are used to weigh objectives in the weight sum method.
                If None, do not multiply BVs by a weight sum vector.

            Returns
            -------
            avgbv : float, numpy.ndarray
                Genotypic mean for the selection.
            """
            # calculate weighted average BVs
            # if traitwt is None:
            #     (j,n).(n,t) -> (j,t)
            # else:
            #     (j,n).{(n,t).(t,) -> (n,)} -> (j,)
            avgbv = sel.dot(bv) if traitwt is None else sel.dot(bv.dot(traitwt))
            return avgbv

        return objfn_vec
