import cvxpy
import math
import numpy
import warnings

from . import SelectionProtocol

import pybropt.core.random
from pybropt.core.error import check_is_int
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_is_Generator

class OptimalContributionSelection(SelectionProtocol):
    """docstring for OptimalContributionSelection."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nparent, ncross, nprogeny, inbfn, cmatcls, bvtype = "gebv", rng = None, **kwargs):
        """
        cmatcls : CoancestryMatrix class
        inbfn : function
            inbfn(t_cur, t_max)
            Returns constraint for mean population inbreeding defined as:
                (1/2) x'Ax = x'Kx <= inbfn(t_cur, t_max)
        bvtype : str
            Whether to use GEBVs or phenotypic EBVs.
            Options:
                "gebv"      Use GEBVs
                "ebv"       Use EBVs
        """
        super(OptimalContributionSelection, self).__init__(**kwargs)

        # error checks
        check_is_int(nparent, "nparent")
        # TODO: check inbfn
        check_is_int(ncross, "ncross")
        check_is_int(nprogeny, "nprogeny")
        # TODO: check cmatcls
        # TODO: check bvtype
        cond_check_is_Generator(rng, "rng")

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.cmatcls = cmatcls
        self.bvtype = bvtype
        self.inbfn = inbfn
        self.rng = pybropt.core.random if rng is None else rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def get_bv(self, pgmat, gmat, bvmat, gpmod, bvtype = None):
        """
        Calculate breeding value matrix for use in optimization.

        Returns
        -------
        out : numpy.ndarray
            Breeding value matrix of shape (n,t).
            Where:
                'n' is the number of individuals.
                't' is the number of traits.
        """
        if bvtype is None:                  # if no bvtype provided
            bvtype = self.bvtype            # get default option
        bvtype = bvtype.lower()             # convert bvtype to lowercase
        if bvtype == "gebv":                # use GEBVs estimated from genomic model
            return gpmod.predict(gmat).mat  # calculate GEBVs
        elif bvtype == "ebv":               # use EBVs estimated by some means
            return bvmat.mat                # get breeding values
        elif bvtype == "tbv":               # use true BVs
            return gpmod.predict(pgmat).mat # calculate true BVs
        else:
            raise ValueError("unknown bvtype: options are 'gebv', 'ebv', 'tbv'")

    def calc_K(self, pgmat, gmat, cmatcls = None, bvtype = None):
        """
        Returns
        -------
        out : numpy.ndarray
            A kinship matrix of shape (n,n).
            Where:
                'n' is the number of individuals.
        """
        # set default parameters
        if cmatcls is None:                     # if no cmatcls provided
            cmatcls = self.cmatcls              # get default
        if bvtype is None:                      # if no bvtype provided
            bvtype = self.bvtype                # get default option
        bvtype = bvtype.lower()                 # convert bvtype to lowercase
        kmat = None                             # decalre kinship matrix variable
        if bvtype == "gebv":                    # use GEBVs estimated from genomic model
            return cmatcls.from_gmat(gmat).mat  # calculate kinship using gmat
        elif bvtype == "ebv":                   # use EBVs estimated by some means
            # TODO: implement pedigree or something
            return cmatcls.from_gmat(gmat).mat  # calculate kinship using gmat
        elif bvtype == "tbv":                   # use true BVs
            return cmatcls.from_gmat(pgmat).mat # calculate true kinship
        else:
            raise ValueError("unknown bvtype: options are 'gebv', 'ebv', 'tbv'")

    def solve_OCS(self, bv, C, inbmax):
        """
        Define and solve OCS using CVXPY.

        Parameters
        ----------
        bv : numpy.ndarray
            Array of shape (n,) containing breeding values for each parent.
        C : numpy.ndarray
            Array of shape (n,n) containing the Cholesky decomposition of the
            kinship matrix. Must be an upper triangle matrix.
        inbmax : float
            Maximum mean inbreeding allowed.

        Returns
        -------
        contrib : numpy.ndarray
            A contribution vector of shape (n,) defining each parent's relative
            contribution.
        """
        # get the number of taxa
        ntaxa = len(bv)

        # define vector variable to optimize
        x = cvxpy.Variable(ntaxa)                   # (n,)

        # define the objective function
        soc_objfn = cvxpy.Maximize(bv @ x)          # max (bv)'(sel)

        # define constraints
        soc_constraints = [
            cvxpy.SOC(inbmax, C @ x),               # ||C @ x||_2 <= inbmax
            cvxpy.sum(x) == 1.0,                    # sum(x_i) == 1
            x >= 0.0                                # x_i >= 0 for all i
        ]

        # define problem
        prob = cvxpy.Problem(
            soc_objfn,                              # maximize yield
            soc_constraints                         # diversity constraint
        )

        # solve the problem
        sol = prob.solve()

        # calculate contributions based on the state of the problem
        contrib = None
        if prob.status != "optimal":                    # if the problem is not optimal, use windowed fitness proportional selection
            contrib = bv - bv.min()                     # window fitnesses
            if contrib.sum() < 1e-10:                   # if everything is near identical
                contrib = numpy.repeat(1/ntaxa, ntaxa)  # default to equal chance
        else:                                           # else, problem has been solved
            contrib = numpy.array(x.value)              # convert solution to numpy.ndarray

        return contrib

    def sus(self, k, contrib):
        """
        Stochastic universal sampling

        k : int
            Number of individuals to sample.
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

    def outcross_shuffle(self, sel):
        """
        Shuffle individuals ensuring they do not mate with themselves.
        """
        self.rng.shuffle(sel)                           # start with random shuffle
        for st in range(0, len(sel), 2):                # for each female
            sp = st+1                                   # get stop index (male)
            j = st+2                                    # get one beyond stop index
            while (sel[st]==sel[sp]) and (j<len(sel)):  # if female == male && within array bounds
                sel[sp], sel[j] = sel[j], sel[sp]       # exchange with next entry
                j += 1                                  # increment beyond index
        return sel

    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, method = "single", nparent = None, ncross = None, nprogeny = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, ndset_trans = None, ndset_trans_kwargs = None, ndset_wt = None, inbfn = None, cmatcls = None, bvtype = None **kwargs):
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        method : str
            Optimization strategy.
            Option  | Description
            --------+-----------------------------------------------------------
            single  | Transform all breeding values into a single overall
                    | breeding value using the function 'objfn_trans'. Then
                    | solve for OCS with a diversity constraint using
                    | transformed breeding values.
            --------+-----------------------------------------------------------
            pareto  | Treat inbreeding and each trait as different objectives.
                    | Transform this list of objectives using 'objfn_trans' to
                    | get a list of transformed objectives. Approximate the
                    | Pareto by identifying a set of non-dominated points along
                    | each transformed objective. Then apply 'ndset_trans' to
                    | score the non-dominated points.
        nparent : int
        ncross : int
        nprogeny : int
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing five objects: (pgmat, sel, ncross, nprogeny, misc)
            pgmat : PhasedGenotypeMatrix
                A PhasedGenotypeMatrix of parental candidates.
            sel : numpy.ndarray
                Array of indices specifying a cross pattern. Each index
                corresponds to an individual in 'pgmat'.
            ncross : numpy.ndarray
                Number of crosses to perform per cross pattern.
            nprogeny : numpy.ndarray
                Number of progeny to generate per cross.
            misc : dict
                Miscellaneous output (user defined).
        """
        # get parameters
        if nparent is None:
            nparent = self.nparent
        if ncross is None:
            ncross = self.ncross
        if nprogeny is None:
            nprogeny = self.nprogeny
        if objfn_trans is None:
            objfn_trans = self.objfn_trans
        if objfn_trans_kwargs is None:
            objfn_trans_kwargs = self.objfn_trans_kwargs
        if objfn_wt is None:
            objfn_wt = self.objfn_wt
        if ndset_trans is None:
            ndset_trans = self.ndset_trans
        if ndset_trans_kwargs is None:
            ndset_trans_kwargs = self.ndset_trans_kwargs
        if ndset_wt is None:
            ndset_wt = self.ndset_wt
        if inbfn is None:
            inbfn = self.inbfn
        if cmatcls is None:
            cmatcls = self.cmatcls
        if bvtype is None:
            bvtype = self.bvtype

        # convert method string to lower
        method = method.lower()

        # Solve OCS using linear programming
        if method == "single":
            ##############################
            # Optimization problem setup #
            ##############################

            # get breeding values
            bv = self.get_bv(pgmat, gmat, bvmat, gpmod)     # (n,t)

            # apply transformation to each breeding value weight
            if objfn_trans:
                # for each row (individual), transform the row to a single objective
                # (n,t) --transform--> (n,)
                bv = numpy.array([objfn_trans(e, **objfn_trans_kwargs) for e in bv])

            # calculate kinship matrix (ndarray)
            K = self.calc_K(pgmat, gmat, cmatcls, bvtype)   # (n,n)

            # get sqrt(max inbreeding)
            inbmax = math.sqrt(inbfn(t_cur, t_max))         # sqrt(max inbreeding)

            # declare contributions variable
            contrib = None

            # cholesky decompose K into K = C'C
            # needed to decompose x'Kx <= inbfn to:
            # ||Cx||_2 <= sqrt(inbfn)
            try:
                C = numpy.linalg.cholesky(K).T              # (n,n)
                contrib = self.solve_OCS(bv, C, inbmax)     # solve OCS
            except numpy.linalg.LinAlgError:
                warnings.warn(
                    "Unable to decompose kinship matrix using Cholesky decomposition: Kinship matrix is not positive definite.\n"+
                    "    This could be caused by lack of genetic diversity.\n"+
                    "Reverting to windowed fitness proportional or uniform selection..."
                )
                ntaxa = len(bv)                             # get number of taxa
                contrib = bv - bv.min()                     # window fitnesses
                if contrib.sum() < 1e-10:                   # if everything is near identical
                    contrib = numpy.repeat(1/ntaxa, ntaxa)  # default to equal chance

            ##################
            # select parents #
            ##################

            # sample selections using stochastic universal sampling
            sel = self.sus(nparent, contrib)                # sample indices

            # make sure parents are forced to outbreed
            sel = self.outcross_shuffle(sel)

            # pack contribution proportions into output dictionary
            misc = {"contrib" : contrib}

            return pgmat, sel, ncross, nprogeny, misc

        # estimate Pareto frontier, then choose from non-dominated points.
        elif method == "pareto":
            # get the pareto frontier
            frontier, sel_config, misc = self.pareto(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                nparent = nparent,
                objfn_trans = objfn_trans,
                objfn_trans_kwargs = objfn_trans_kwargs,
                objfn_wt = objfn_wt
            )

            # get scores for each of the points along the pareto frontier
            score = ndset_wt * ndset_trans(frontier, **ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to misc
            misc["frontier"] = frontier
            misc["sel_config"] = sel_config

            return pgmat, sel_config[ix], ncross, nprogeny, misc

        # raise error otherwise
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, cmatcls = None, bvtype = None, **kwargs):
        """
        Return a parent selection objective function.
        """
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs
        if cmatcls is None:
            cmatcls = self.cmatcls
        if bvtype is None:
            bvtype = self.bvtype

        # get pointers to raw numpy.ndarray matrices
        mat = self.get_bv(pgmat, gmat, bvmat, gpmod)    # (n,t) get breeding values
        K = self.calc_K(pgmat, gmat, cmatcls, bvtype)   # (n,n) get kinship matrix

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,     # byte code pointer
            self.objfn_static.__globals__,  # global variables
            None,                           # new name for the function
            (mat, K, trans, trans_kwargs),  # default values for arguments
            self.objfn_static.__closure__   # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, trans = None, trans_kwargs = None, cmatcls = None, bvtype = None, **kwargs):
        """
        Return a vectorized objective function.
        """
        # get default parameters if any are None
        if trans is None:
            trans = self.objfn_trans
        if trans_kwargs is None:
            trans_kwargs = self.objfn_trans_kwargs
        if cmatcls is None:
            cmatcls = self.cmatcls
        if bvtype is None:
            bvtype = self.bvtype

        # get pointers to raw numpy.ndarray matrices
        mat = self.get_bv(pgmat, gmat, bvmat, gpmod)    # (n,t) get breeding values
        K = self.calc_K(pgmat, gmat, cmatcls, bvtype)   # (n,n) get kinship matrix

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (mat, K, trans, trans_kwargs),      # default values for arguments
            self.objfn_vec_static.__closure__   # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, nparent = None, objfn_trans = None, objfn_trans_kwargs = None, objfn_wt = None, **kwargs):
        """
        Calculate a Pareto frontier for objectives.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotypes
        ptdf : PhenotypeDataFrame
            Phenotype dataframe
        bvmat : BreedingValueMatrix
            Breeding value matrix
        gpmod : GenomicModel
            Genomic prediction model
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing three objects (frontier, sel_config, misc)
            Elements
            --------
            frontier : numpy.ndarray
                Array of shape (q,v) containing Pareto frontier points.
                Where:
                    'q' is the number of points in the frontier.
                    'v' is the number of objectives for the frontier.
            sel_config : numpy.ndarray
                Array of shape (q,k) containing parent selection decisions for
                each corresponding point in the Pareto frontier.
                Where:
                    'q' is the number of points in the frontier.
                    'k' is the number of search space decision variables.
            misc : dict
                A dictionary of miscellaneous output. (User specified)
        """
        raise NotImplementedError("method is abstract")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, K, trans, kwargs):
        """
        Score a parent contribution vector according to its expected breeding
        value.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A parent contribution vector of shape (n,) and floating dtype.
            Where:
                'n' is the number of individuals.
        mat : numpy.ndarray
            A breeding value matrix of shape (n,t).
            Where:
                'n' is the number of individuals.
                't' is the number of traits.
        K : numpy.ndarray
            A kinship matrix of shape (n,n).
            Where:
                'n' is the number of individuals.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        ocs : numpy.ndarray
            A EBV matrix of shape (1+t,) if 'trans' is None.
            The first index in the array is the mean expected kinship:
                mean expected inbreeding = (sel') @ K @ (sel)
            Other indices are the mean expected trait values for the other 't'
            traits.
            Otherwise, of shape specified by 'trans'.
            Where:
                't' is the number of traits.
        """
        # Calculate the mean expected kinship: x'Kx
        # Step 1: (n,) . (n,n) -> (n,)
        # Step 2: (n,) . (n,) -> scalar
        inb = sel.dot(K).dot(sel)

        # OCS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # calculate mean expected BVs
        mebv = sel.dot(mat)

        # append values together with inbreeding first
        # scalar and (t,) --append--> (1+t,)
        ocs = numpy.insert(mebv, 0, inb)

        # apply transformations
        # (1+t,) ---trans---> (?,)
        if trans:
            ocs = trans(ocs, **kwargs)

        return ocs

    @staticmethod
    def objfn_vec_static(sel, mat, K, trans, kwargs):
        """
        Score a parent contribution vector according to its expected breeding
        value.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A parent contribution vector of shape (j,n) and floating dtype.
            Where:
                'j' is the number of selection configurations.
                'n' is the number of individuals.
        mat : numpy.ndarray
            A breeding value matrix of shape (n,t).
            Where:
                'n' is the number of individuals.
                't' is the number of traits.
        K : numpy.ndarray
            A kinship matrix of shape (n,n).
            Where:
                'n' is the number of individuals.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:
                Must accept a single numpy.ndarray argument.
                Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to 'trans' function.

        Returns
        -------
        cgs : numpy.ndarray
            A EBV matrix of shape (j,1+t) if 'trans' is None.
            The first column in the matrix is the mean expected kinship:
                mean expected inbreeding = (sel') @ K @ (sel)
            Other indices are the mean expected trait values for the other 't'
            traits.
            Otherwise, of shape specified by 'trans'.
            Where:
                'j' is the number of selection configurations.
                't' is the number of traits.
        """
        # Calculate the mean expected kinship: x'Kx
        # Step 1: for each row in range {1,2,...,j}
        # Step 2: (n,) . (n,n) -> (n,)
        # Step 3: (n,) . (n,) -> scalar
        # Step 4: (j,)                              # construct array
        inb = numpy.array([e.dot(K).dot(e) for e in sel])

        # OCS calculation explanation
        # Step 1: (j,n) @ (n,t) -> (j,t)    # calculate mean expected BVs
        mebv = sel @ mat

        # append values together with inbreeding first
        # (j,) and (j,t) --append--> (j,1+t)
        ocs = numpy.insert(mebv, 0, inb, axis = 1)

        # apply transformations
        # (j,1+t) ---trans---> (?,?)
        if trans:
            ocs = trans(ocs, **kwargs)

        return ocs
