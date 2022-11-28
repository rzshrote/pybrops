"""
Module implementing selection protocols for optimal contribution selection.
"""

from distutils.log import warn
from optparse import Option
import cvxpy
import math
import numpy
import warnings
import types
from typing import Callable
from typing import Union
from typing import Optional

from pybrops.algo.opt.NSGA3UnityConstraintGeneticAlgorithm import NSGA3UnityConstraintGeneticAlgorithm
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_inherits
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.random import global_prng
from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.cmat.DenseVanRadenCoancestryMatrix import DenseVanRadenCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.breed.prot.sel.sampling import stochastic_universal_sampling
from pybrops.breed.prot.sel.sampling import two_way_outcross_shuffle

class ContinuousOptimalContributionSelection(SelectionProtocol):
    """
    Class implementing selection protocols for genomic optimal contribution selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nparent: int, ncross: int, nprogeny: int, inbfn: Callable,
        cmatcls: type = DenseVanRadenCoancestryMatrix, bvtype = "gebv", method = "single",
        objfn_trans: Callable = None, objfn_trans_kwargs: dict = None, objfn_wt: Union[float,numpy.ndarray] = 1.0,
        ndset_trans: Callable = None, ndset_trans_kwargs: dict = None, ndset_wt: float = -1.0,
        moalgo = None, rng = global_prng, **kwargs):
        """
        Constructor for Optimal Contribution Selection (OCS).

        Parameters
        ----------
        nparent : int
            Number of parents to select.
        ncross : int
            Number of crosses per configuration.
        nprogeny : int
            Number of progeny to derive from each cross.
        inbfn : function
            Inbreeding control function: ``inbfn(t_cur, t_max)``.

            Returns constraint for mean population inbreeding defined as:

            .. math::

                \\frac{1}{2} \\textbf{x}' \\textbf{A} \\textbf{x} =
                \\textbf{x}' \\textbf{K} \\textbf{x} =
                f^{\\textup{Inb}}(t_{cur}, t_{max})

            Where:

            - :math:`x` is the parental contribution vector.
            - :math:`A` is the additive relationship matrix.
            - :math:`K` is the kinship relationship matrix.
            - :math:`f^{\\textup{Inb}}` is ``inbfn``.
            - :math:`t_{cur}` is the current time.
            - :math:`t_{max}` is the deadline time.
        cmatcls : class
            The class name of a CoancestryMatrix to generate.
        bvtype : str
            Whether to use GEBVs or phenotypic EBVs.

            +------------+-------------+
            | Option     | Description |
            +============+=============+
            | ``"gebv"`` | Use GEBVs   |
            +------------+-------------+
            | ``"ebv"``  | Use EBVs    |
            +------------+-------------+
        method : str
            Optimization strategy.

            +--------+---------------------------------------------------------+
            | Option | Description                                             |
            +========+=========================================================+
            | single | Transform all breeding values into a single overall     |
            |        | breeding value using the function ``objfn_trans``. Then |
            |        | solve for OCS with a diversity constraint using         |
            |        | transformed breeding values.                            |
            +--------+---------------------------------------------------------+
            | pareto | Treat inbreeding and each trait as different            |
            |        | objectives. Transform this list of objectives using     |
            |        | ``objfn_trans`` to get a list of transformed            |
            |        | objectives. Approximate the Pareto by identifying a set |
            |        | of non-dominated points along each transformed          |
            |        | objective. Then apply ``ndset_trans`` to score the      |
            |        | non-dominated points.                                   |
            +--------+---------------------------------------------------------+
        objfn_trans : function or callable
            Function to transform the OCS objective function.

            If method = "single", this function must accept an array of length
            ``(t,)`` and return a scalar, where ``t`` is the number of trait
            breeding values for an individual.

            If method = "pareto", this function must accept an array of length
            ``(1+t,)`` and return a numpy.ndarray, where ``t`` is the number of
            trait breeding values for an individual. The first element of the
            input array is the mean inbreeding coefficient for the selection.

            General function definition::

                objfn_trans(obj, **kwargs):
                    Parameters
                        obj : scalar, numpy.ndarray
                            Objective scalar or vector to be transformed
                        kwargs : dict
                            Additional keyword arguments
                    Returns
                        out : scalar, numpy.ndarray
                            Transformed objective scalar or vector.
        objfn_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'objfn_trans'.
        objfn_wt : float, numpy.ndarray
            Weight applied to transformed objective function. Indicates whether
            a function is maximizing or minimizing:

            - ``1.0`` for maximizing function.
            - ``-1.0`` for minimizing function.
        ndset_trans : numpy.ndarray
            Function to transform nondominated points along the Pareto frontier
            into a single score for each point.

            Function definition::

                ndset_trans(ndset, **kwargs):
                    Parameters
                        ndset : numpy.ndarray
                            Array of shape (j,o) containing nondominated points.
                            Where 'j' is the number of nondominated points and
                            'o' is the number of objectives.
                        kwargs : dict
                            Additional keyword arguments.
                    Returns
                        out : numpy.ndarray
                            Array of shape (j,) containing transformed Pareto
                            frontier points.
        ndset_trans_kwargs : dict
            Dictionary of keyword arguments to be passed to 'ndset_trans'.
        ndset_wt : float
            Weight applied to transformed nondominated points along Pareto
            frontier. Indicates whether a function is maximizing or minimizing.
                1.0 for maximizing function.
                -1.0 for minimizing function.
        moalgo : OptimizationAlgorithm
            Multi-objective optimization algorithm to optimize the objective
            functions. If ``None``, use a NSGA3UnityConstraintGeneticAlgorithm
            with the following parameters::

                moalgo = NSGA3UnityConstraintGeneticAlgorithm(
                    ngen = 250,             # number of generations to evolve
                    mu = 100,               # number of parents in population
                    lamb = 100,             # number of progeny to produce
                    cxeta = 30.0,           # crossover variance parameter
                    muteta = 20.0,          # mutation crossover parameter
                    refpnts = None,         # hyperplane reference points
                    save_logbook = False,   # whether to save logs or not
                    rng = self.rng          # PRNG source
                )
        rng : numpy.random.Generator or None
            A random number generator source. Used for optimization algorithms.
        """
        super(ContinuousOptimalContributionSelection, self).__init__(**kwargs)

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.inbfn = inbfn
        self.cmatcls = cmatcls
        self.bvtype = bvtype
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = rng
        self.moalgo = moalgo    # must go after rng initialization!!!

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def nparent():
        doc = "The nparent property."
        def fget(self):
            return self._nparent
        def fset(self, value):
            check_is_int(value, "nparent")      # must be int
            check_is_gt(value, "nparent", 0)    # int must be >0
            self._nparent = value
        def fdel(self):
            del self._nparent
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    nparent = property(**nparent())

    def ncross():
        doc = "The ncross property."
        def fget(self):
            return self._ncross
        def fset(self, value):
            check_is_int(value, "ncross")       # must be int
            check_is_gt(value, "ncross", 0)     # int must be >0
            self._ncross = value
        def fdel(self):
            del self._ncross
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ncross = property(**ncross())

    def nprogeny():
        doc = "The nprogeny property."
        def fget(self):
            return self._nprogeny
        def fset(self, value):
            check_is_int(value, "nprogeny")     # must be int
            check_is_gt(value, "nprogeny", 0)   # int must be >0
            self._nprogeny = value
        def fdel(self):
            del self._nprogeny
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    nprogeny = property(**nprogeny())

    def inbfn():
        doc = "Inbreeding control function."
        def fget(self):
            return self._inbfn
        def fset(self, value):
            check_is_callable(value, "inbfn")
            self._inbfn = value
        def fdel(self):
            del self._inbfn
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    inbfn = property(**inbfn())

    def cmatcls():
        doc = "Coancestry matrix class."
        def fget(self):
            return self._cmatcls
        def fset(self, value):
            check_inherits(value, "cmatcls", CoancestryMatrix)
            self._cmatcls = value
        def fdel(self):
            del self._cmatcls
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    cmatcls = property(**cmatcls())

    def bvtype():
        doc = "Breeding value matrix type."
        def fget(self):
            return self._bvtype
        def fset(self, value):
            check_is_str(value, "bvtype")   # must be string
            value = value.lower()           # convert to lowercase
            options = ("gebv", "ebv")       # method options
            if value not in options:            # if not method supported
                raise ValueError(               # raise ValueError
                    "Unsupported 'method'. Options are: " +
                    ", ".join(map(str, options))
                )
            self._bvtype = value
        def fdel(self):
            del self._bvtype
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    bvtype = property(**bvtype())

    def method():
        doc = "The method property."
        def fget(self):
            return self._method
        def fset(self, value):
            check_is_str(value, "method")       # must be string
            value = value.lower()               # convert to lowercase
            options = ("single", "pareto")      # method options
            if value not in options:            # if not method supported
                raise ValueError(               # raise ValueError
                    "Unsupported 'method'. Options are: " +
                    ", ".join(map(str, options))
                )
            self._method = value
        def fdel(self):
            del self._method
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    method = property(**method())

    def objfn_trans():
        doc = "The objfn_trans property."
        def fget(self):
            return self._objfn_trans
        def fset(self, value):
            if value is not None:                       # if given object
                check_is_callable(value, "objfn_trans") # must be callable
            self._objfn_trans = value
        def fdel(self):
            del self._objfn_trans
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    objfn_trans = property(**objfn_trans())

    def objfn_trans_kwargs():
        doc = "The objfn_trans_kwargs property."
        def fget(self):
            return self._objfn_trans_kwargs
        def fset(self, value):
            if value is None:                           # if given None
                value = {}                              # set default to empty dict
            check_is_dict(value, "objfn_trans_kwargs")  # check is dict
            self._objfn_trans_kwargs = value
        def fdel(self):
            del self._objfn_trans_kwargs
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    objfn_trans_kwargs = property(**objfn_trans_kwargs())

    def objfn_wt():
        doc = "The objfn_wt property."
        def fget(self):
            return self._objfn_wt
        def fset(self, value):
            self._objfn_wt = value
        def fdel(self):
            del self._objfn_wt
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    objfn_wt = property(**objfn_wt())

    def ndset_trans():
        doc = "The ndset_trans property."
        def fget(self):
            return self._ndset_trans
        def fset(self, value):
            if value is not None:                       # if given object
                check_is_callable(value, "ndset_trans") # must be callable
            self._ndset_trans = value
        def fdel(self):
            del self._ndset_trans
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ndset_trans = property(**ndset_trans())

    def ndset_trans_kwargs():
        doc = "The ndset_trans_kwargs property."
        def fget(self):
            return self._ndset_trans_kwargs
        def fset(self, value):
            if value is None:                           # if given None
                value = {}                              # set default to empty dict
            check_is_dict(value, "ndset_trans_kwargs")  # check is dict
            self._ndset_trans_kwargs = value
        def fdel(self):
            del self._ndset_trans_kwargs
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ndset_trans_kwargs = property(**ndset_trans_kwargs())

    def ndset_wt():
        doc = "The ndset_wt property."
        def fget(self):
            return self._ndset_wt
        def fset(self, value):
            self._ndset_wt = value
        def fdel(self):
            del self._ndset_wt
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    ndset_wt = property(**ndset_wt())

    def moalgo():
        doc = "The moalgo property."
        def fget(self):
            return self._moalgo
        def fset(self, value):
            if value is None:
                value = NSGA3UnityConstraintGeneticAlgorithm(
                    ngen = 600,             # number of generations to evolve
                    mu = 100,               # number of parents in population
                    lamb = 100,             # number of progeny to produce
                    cxeta = 30.0,           # crossover variance parameter
                    muteta = 20.0,          # mutation crossover parameter
                    refpnts = None,         # hyperplane reference points
                    save_logbook = False,   # whether to save logs or not
                    rng = self.rng          # PRNG source
                )
            self._moalgo = value
        def fdel(self):
            del self._moalgo
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    moalgo = property(**moalgo())

    def rng():
        doc = "The rng property."
        def fget(self):
            return self._rng
        def fset(self, value):
            if value is None:               # if rng is None
                value = global_prng         # use default random number generator
            check_is_Generator_or_RandomState(value, "rng")# check is numpy.Generator
            self._rng = value
        def fdel(self):
            del self._rng
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    rng = property(**rng())

    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
    def _get_bv(
            self, 
            pgmat: Optional[PhasedGenotypeMatrix], 
            gmat: Optional[GenotypeMatrix], 
            bvmat: Optional[BreedingValueMatrix], 
            gpmod: Optional[GenomicModel]
        ) -> BreedingValueMatrix:
        """
        Calculate breeding value matrix for use in optimization.

        Returns
        -------
        out : numpy.ndarray
            Breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        """
        if self.bvtype == "gebv":       # use GEBVs estimated from genomic model
            return gpmod.gebv(gmat)     # calculate GEBVs
        elif self.bvtype == "ebv":      # use EBVs estimated by some means
            return bvmat                # get breeding values
        elif self.bvtype == "tbv":      # use true BVs
            return gpmod.gebv(pgmat)    # calculate true BVs

    def _calc_G(
            self, 
            pgmat: Optional[PhasedGenotypeMatrix], 
            gmat: Optional[GenotypeMatrix]
        ) -> CoancestryMatrix:
        """
        Calculate a kinship matrix from a genotype matrix.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix, None
            True phased genome matrix.
        gmat : GenotypeMatrix, None
            Genotype matrix.

        Returns
        -------
        out : CoancestryMatrix
            A coancestry matrix of shape ``(n,n)``.

            Where:

            - ``n`` is the number of individuals.
        """
        if self.bvtype == "gebv":                   # use GEBVs estimated from genomic model
            return self.cmatcls.from_gmat(gmat)     # calculate kinship using gmat
        elif self.bvtype == "ebv":                  # use EBVs estimated by some means
            # TODO: implement pedigree or something
            return self.cmatcls.from_gmat(gmat)     # calculate kinship using gmat
        elif self.bvtype == "tbv":                  # use true BVs
            return self.cmatcls.from_gmat(pgmat)    # calculate true kinship

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs):
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
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing four objects: ``(pgmat, sel, ncross, nprogeny)``.

            Where:

            - ``pgmat`` is a PhasedGenotypeMatrix of parental candidates.
            - ``sel`` is a ``numpy.ndarray`` of indices specifying a cross
              pattern. Each index corresponds to an individual in ``pgmat``.
            - ``ncross`` is a ``numpy.ndarray`` specifying the number of
              crosses to perform per cross pattern.
            - ``nprogeny`` is a ``numpy.ndarray`` specifying the number of
              progeny to generate per cross.
        """
        # Solve OCS using second order cone programming
        if self.method == "single":
            ##############################
            # Optimization problem setup #
            ##############################

            # get breeding values: (n,t)
            bv = self._get_bv(pgmat, gmat, bvmat, gpmod).mat

            # apply transformation to each breeding value weight
            if self.objfn_trans:
                # for each row (individual), transform the row to a single objective
                # (n,t) --transform--> (n,)
                bv = numpy.array([self.objfn_trans(e, **self.objfn_trans_kwargs) for e in bv])

            # make sure transformation is done correctly.
            if bv.ndim > 1:
                raise RuntimeError("objfn_trans does not return a scalar value")

            # get genomic relationship matrix: (n,n)
            G = self._calc_G(pgmat, gmat)

            # declare contributions variable
            contrib = None

            # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
            # if we successfully were able to apply the jitter, then perform optimization.
            if G.apply_jitter():
                K = G.mat_asformat("kinship")                   # convert G to (1/2)G (kinship analogue): (n,n)
                K_inv = numpy.linalg.inv(K)                     # get K^-1 needed for bounds check

                inb_target = self.inbfn(t_cur, t_max)           # get target inbreeding level
                inb_min = 1.0 / K_inv.sum()                     # 1 / (1'K^(-1)1)
                inb_max = numpy.max(numpy.trace(K))             # max(trace(K))

                # test for infeasibility
                if inb_target < inb_min:
                    warnings.warn(
                        "Provided inbreeding target of {0} is infeasible.\n".format(inb_target) +
                        "Increasing inbreeding target to {0} which is in feasible region...".format(inb_min + 1e-6)
                    )
                    # give some wiggle room for numerical inaccuracies due to matrix
                    # inversion, which is by nature unstable.
                    inb_target = inb_min + 1e-6
                if inb_target > inb_max:
                    warnings.warn(
                        "Provided inbreeding target of {0} is infeasible.\n".format(inb_target) +
                        "Decreasing inbreeding target to {0} which is in feasible region...".format(inb_max)
                    )
                    inb_target = inb_max

                # calculate constants for optimization
                C = numpy.linalg.cholesky(K).T                  # cholesky decomposition of K matrix: (n,n)

                # get the number of taxa
                ntaxa = len(bv)

                # define vector variable to optimize
                solution = cvxpy.Variable(ntaxa)                    # (n,)

                # define the objective function
                soc_objfn = cvxpy.Maximize(bv @ solution)           # max (bv)'(sel)

                # define constraints
                soc_constraints = [
                    cvxpy.SOC(math.sqrt(inb_target), C @ solution), # ||C @ x||_2 <= sqrt(max inbreeding)
                    cvxpy.sum(solution) == 1.0,                     # sum(x_i) == 1
                    solution >= 0.0                                 # x_i >= 0 for all i
                ]

                # define problem
                problem = cvxpy.Problem(
                    soc_objfn,                              # maximize yield
                    soc_constraints                         # diversity constraint
                )

                # solve the problem
                problem.solve()

                # store solution results
                if problem.status == "optimal":
                    contrib = numpy.array(solution.value)       # convert solution to numpy.ndarray
                else:
                    warnings.warn("Problem.status == {0}\n".format(problem.status))
            else:
                warnings.warn(
                    "Unable to solve SOCP: Kinship matrix is not positive definite.\n"+
                    "    This could be caused by lack of genetic diversity.\n"
                )
            
            if contrib is None:
                warnings.warn("Reverting to windowed fitness proportional or uniform selection...")
                ntaxa = len(bv)                                 # get number of taxa
                contrib = bv - bv.min()                         # window fitnesses
                if contrib.sum() < 1e-10:                       # if everything is near identical
                    contrib = numpy.repeat(1/ntaxa, ntaxa)      # default to equal chance

            ##################
            # select parents #
            ##################

            # sample selections using stochastic universal sampling
            sel = stochastic_universal_sampling(self.nparent, contrib, self.rng)

            # make sure parents are forced to outbreed
            sel = two_way_outcross_shuffle(sel, self._rng)

            # pack contribution proportions into output dictionary
            if miscout is not None:
                miscout["contrib"] = contrib

            return pgmat, sel, self.ncross, self.nprogeny

        # estimate Pareto frontier, then choose from non-dominated points.
        elif self.method == "pareto":
            # get the pareto frontier
            frontier, sel_config = self.pareto(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                miscout = miscout,
                **kwargs
            )

            # get scores for each of the points along the pareto frontier
            score = self.ndset_wt * self.ndset_trans(frontier, **self.ndset_trans_kwargs)

            # get index of maximum score
            ix = score.argmax()

            # add fields to miscout
            if miscout is not None:
                miscout["frontier"] = frontier
                miscout["sel_config"] = sel_config

            return pgmat, sel_config[ix], self.ncross, self.nprogeny

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        """
        Return an objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Input phased genotype matrix containing genomes.
        gmat : GenotypeMatrix
            Input genotype matrix.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Input breeding value matrix.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A selection objective function for the specified problem.
        """
        # get default parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # get pointers to raw numpy.ndarray matrices
        mat = self._get_bv(pgmat, gmat, bvmat, gpmod).mat   # (n,t) get breeding values
        G = self._calc_G(pgmat, gmat)                       # (n,n) get genomic relationship matrix
        K = G.mat_asformat("kinship")                       # (n,n) kinship matrix

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

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs):
        """
        Return a vectorized objective function.
        """
        # get default parameters if any are None
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # get pointers to raw numpy.ndarray matrices
        mat = self._get_bv(pgmat, gmat, bvmat, gpmod).mat   # (n,t) get breeding values
        G = self._calc_G(pgmat, gmat)                       # (n,n) get genomic relationship matrix
        K = G.mat_asformat("kinship")                       # (n,n) kinship matrix

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

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs):
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing two objects ``(frontier, sel_config)``.

            Where:

            - ``frontier`` is a ``numpy.ndarray`` of shape ``(q,v)`` containing
              Pareto frontier points.
            - ``sel_config`` is a ``numpy.ndarray`` of shape ``(q,k)`` containing
              parent selection decisions for each corresponding point in the
              Pareto frontier.

            Where:

            - ``q`` is the number of points in the frontier.
            - ``v`` is the number of objectives for the frontier.
            - ``k`` is the number of search space decision variables.
        """
        # get selection parameters
        objfn_wt = self.objfn_wt

        # get number of taxa
        ntaxa = gmat.ntaxa

        # create objective function
        objfn = self.objfn(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max
        )

        # create search space
        sspace = numpy.stack(
            [numpy.repeat(0.0, ntaxa), numpy.repeat(1.0, ntaxa)]
        )

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,          # objective function
            k = ntaxa,              # vector length to optimize (sspace^k)
            sspace = sspace,        # search space options
            objfn_wt = objfn_wt,    # weights to apply to each objective
            **kwargs
        )

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        return frontier, sel_config

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
            A parent contribution vector of shape ``(n,)`` and floating dtype.

            Where:

            - ``n`` is the number of individuals.
        mat : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        K : numpy.ndarray
            A kinship matrix of shape (n,n).

            Where:

            - ``n`` is the number of individuals.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        ocs : numpy.ndarray
            A EBV matrix of shape (1+t,) if ``trans`` is ``None``.

            The first index in the array is the mean expected kinship:

            .. math::
                mean expected inbreeding = \\textbf{(sel)}' \\textbf{K} \\textbf{(sel)}

            Other indices are the mean expected trait values for the other ``t``
            traits. Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
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
            A parent contribution vector of shape ``(j,n)`` and floating dtype.

            Where:

            - ``j`` is the number of selection configurations.
            - ``n`` is the number of individuals.
        mat : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        K : numpy.ndarray
            A kinship matrix of shape ``(n,n)``.

            Where:

            - ``n`` is the number of individuals.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single numpy.ndarray argument.
            - Must return a single object, whether scalar or numpy.ndarray.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cgs : numpy.ndarray
            A EBV matrix of shape ``(j,1+t)`` if ``trans`` is ``None``.
            The first column in the matrix is the mean expected kinship:

            .. math::
                mean expected inbreeding = \\textbf{(sel)}' \\textbf{K} \\textbf{(sel)}

            Other indices are the mean expected trait values for the other ``t``
            traits. Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.
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
