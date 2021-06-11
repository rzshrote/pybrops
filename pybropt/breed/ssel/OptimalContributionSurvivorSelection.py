import cvxpy
import math

from . import SurvivorSelectionOperator

import pybropt.core.random
from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import cond_check_is_Generator

class OptimalContributionSurvivorSelection(SurvivorSelectionOperator):
    """docstring for OptimalContributionSurvivorSelection."""

    def __init__(self, k_s, traitwt_s, inbfn_s, cmatcls, bvtype = "gebv", rng = None, **kwargs):
        super(OptimalContributionSurvivorSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(k_s, "k_s")
        cond_check_is_ndarray(traitwt_s, "traitwt_s")
        # TODO: check inbfn_s
        # TODO: check cmatcls
        # TODO: check bvtype
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.k_s = k_s
        self.traitwt_s = traitwt_s
        self.cmatcls = cmatcls
        self.bvtype = bvtype
        self.rng = pybropt.core.random if rng is None else rng
        raise NotImplementedError("under construction")

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def get_bv(self, geno, bval, gmod, **kwargs):
        option = self.bvtype.lower()
        if option == "gebv":
            return gmod["main"].predict(geno["main"]).mat   # (n,t)
        elif option == "ebv":
            return bval["main"].mat                         # (n,t)
        else:
            raise RuntimeError("unknown bvtype")

    def calc_K(self, gmat, **kwargs):
        kmat = cmatcls.from_gmat(gmat)  # get coancestry (kinship) matrix
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
        offset = random.uniform(0.,ptr_dist)            # get random start point
        sel = []                                        # declare output list
        ix = 0                                          # index for cumsum
        ptrs = numpy.arange(offset, tot_fit, ptr_dist)  # create pointers
        for ptr in ptrs:                                # for each pointer
            while cumsum[ix] < ptr:                     # advance index to correct location
                ix += 1
            sel.append(indices[ix])                     # append index with element
        sel = numpy.array(sel)                          # convert to ndarray
        return sel

    def sselect(self, t_cur, t_max, geno, bval, gmod, k = None, traitwt = None, **kwargs):
        """
        Select survivors to serve as potential parents for breeding.

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
        trait : numpy.ndarray, None

        Returns
        -------
        out : tuple
            A tuple containing four objects: (geno_new, bval_new, gmod_new, misc)
            geno_new : dict
                A dict containing genotypic data for all breeding populations.
                Must have the following fields:
                    Field | Type                         | Description
                    ------+------------------------------+----------------------
                    cand  | PhasedGenotypeMatrix         | Parental candidate breeding population
                    main  | PhasedGenotypeMatrix         | Main breeding population
            bval_new : dict
                A dict containing breeding value data.
                Must have the following fields:
                    Field      | Type                        | Description
                    -----------+-----------------------------+------------------
                    cand       | BreedingValueMatrix         | Parental candidate breeding population breeding values
                    cand_true  | BreedingValueMatrix         | Parental candidate population true breeding values
                    main       | BreedingValueMatrix         | Main breeding population breeding values
                    main_true  | BreedingValueMatrix         | Main breeding population true breeding values
            gmod_new : dict
                A dict containing genomic models.
                Must have the following fields:
                    Field | Type                 | Description
                    ------+----------------------+------------------------------
                    cand  | GenomicModel         | Parental candidate breeding population genomic model
                    main  | GenomicModel         | Main breeding population genomic model
                    true  | GenomicModel         | True genomic model for trait(s)
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if k is None:
            k = self.k_s
        if traitwt is None:
            traitwt = self.traitwt_s

        # TODO: implement multiobjective/non-dominated OCS
        if traitwt is None:
            raise RuntimeError("multi-objective optimal contribution selection not implemented")

        ##############################
        # Optimization problem setup #
        ##############################
        # TODO: probably put this optimization into a single function, or a
        # static library for code reusability
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

        # do not need to shuffle selections...
        # shallow copy dict
        geno_new = dict(geno)
        bval_new = dict(bval)
        gmod_new = dict(gmod)

        # update cand fields
        geno_new["cand"] = geno["main"].select(sel, axis = 1)
        bval_new["cand"] = bval["main"].select(sel, axis = 0)
        bval_new["cand_true"] = bval["main_true"].select(sel, axis = 0)
        gmod_new["cand"] = gmod["main"]

        misc = {}   # empty dictionary

        return geno_new, bval_new, gmod_new, misc

    def sobjfn(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
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

    def sobjfn_vec(self, t_cur, t_max, geno, bval, gmod, traitwt = None, **kwargs):
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
