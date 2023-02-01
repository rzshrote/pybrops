"""
Module implementing selection protocols for minimum mean genomic relationship selection.
"""

import cvxpy
import numpy
import warnings
import types
from typing import Optional
from typing import Callable
from typing import Type

from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.error import check_inherits
from pybrops.core.error import check_is_class
from pybrops.core.random import global_prng
from pybrops.popgen.cmat.CoancestryMatrix import CoancestryMatrix
from pybrops.popgen.cmat.DenseVanRadenCoancestryMatrix import DenseVanRadenCoancestryMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.breed.prot.sel.sampling import stochastic_universal_sampling
from pybrops.breed.prot.sel.sampling import two_way_outcross_shuffle
from pybrops.algo.opt.MemeticSetGeneticAlgorithm import MemeticSetGeneticAlgorithm

class BinaryMinimumMeanGenomicRelationshipSelection(SelectionProtocol):
    """
    Class implementing selection protocols for minimum mean genomic relationship selection.

    Minimum mean genomic relationship selection (MMGRS) is defined as:

    .. math::
        \\min_{\\mathbf{x}} f_{MMGRS}(\\mathbf{x}) = \\frac{1}{2}\\mathbf{x'Gx}

    With constraints:

    .. math::
        \\mathbf{1_{n}'x} = 1

        \\mathbf{x} \\in \\mathbb{R}^n

    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, 
        nparent: int, 
        ncross: int, 
        nprogeny: int,
        gtype: str = "gmat", 
        gcls: Type[CoancestryMatrix] = DenseVanRadenCoancestryMatrix, 
        method: str = "single",
        objfn_trans: Callable = None, 
        objfn_trans_kwargs: dict = None, 
        objfn_wt = -1.0,
        rng = global_prng, 
        soalgo = None, 
        **kwargs: dict
        ):
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
        gtype : str
            Whether to use genotypes or genomes to calculate genomic relationship.

            +-------------+---------------+
            | Option      | Description   |
            +=============+===============+
            | ``"gmat"``  | Use genotypes |
            +-------------+---------------+
            | ``"pgmat"`` | Use genomes   |
            +-------------+---------------+
        gcls : Any
            Which genomic relationship class to use to construct the genomic relationship matrix.
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

                objfn_trans(obj, **kwargs: dict):
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

                ndset_trans(ndset, **kwargs: dict):
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
        super(BinaryMinimumMeanGenomicRelationshipSelection, self).__init__(**kwargs)
        
        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.gtype = gtype
        self.gcls = gcls
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.rng = rng
        self.soalgo = soalgo

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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    nprogeny = property(**nprogeny())

    def gtype():
        doc = "Genomic relationship matrix type"
        def fget(self):
            return self._gtype
        def fset(self, value):
            check_is_str(value, "gtype")
            value = value.lower()
            if value not in ["gmat","pgmat"]:
                raise ValueError("'gtype' must be 'gmat' or 'pgmat'")
            self._gtype = value
        def fdel(self):
            del self._gtype
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    gtype = property(**gtype())

    def gcls():
        doc = "Genomic relationship matrix class"
        def fget(self):
            return self._gcls
        def fset(self, value):
            check_is_class(value, "gcls")
            check_inherits(value, "gcls", CoancestryMatrix)
            self._gcls = value
        def fdel(self):
            del self._gcls
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    gcls = property(**gcls())

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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    objfn_trans_kwargs = property(**objfn_trans_kwargs())

    def objfn_wt():
        doc = "The objfn_wt property."
        def fget(self):
            return self._objfn_wt
        def fset(self, value):
            self._objfn_wt = value
        def fdel(self):
            del self._objfn_wt
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    objfn_wt = property(**objfn_wt())

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
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    rng = property(**rng())

    def soalgo():
        doc = "The soalgo property."
        def fget(self):
            """Get value for soalgo."""
            return self._soalgo
        def fset(self, value):
            """Set value for soalgo."""
            if value is None:
                value = MemeticSetGeneticAlgorithm(
                    ngen = 100,
                    mu = 100,
                    lamb = 100,
                    M = 1.5,
                    rng = self.rng
                )
            self._soalgo = value
        def fdel(self):
            """Delete value for soalgo."""
            del self._soalgo
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    soalgo = property(**soalgo())

    ############################################################################
    ########################## Private Object Methods ##########################
    ############################################################################
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
        # use genotypes to determine kinship
        if self.gtype == "gmat":
            return self.gcls.from_gmat(gmat)
        # use genomes to determine kinshipe
        elif self.gtype == "pgmat":
            return self.gcls.from_gmat(pgmat)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
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
        # Solve problem using quadratic programming
        if self.method == "single":
            # get vectorized objective function
            objfn = self.objfn(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max
            )

            # optimize using single objective algorithm
            sel_score, sel, misc = self.soalgo.optimize(
                objfn,                              # objective function
                k = self.nparent,                   # number of parents to select
                sspace = numpy.arange(pgmat.ntaxa), # parental indices
                objfn_wt = self.objfn_wt,           # maximizing function
                **kwargs
            )

            # shuffle selection to ensure random mating
            numpy.random.shuffle(sel)

            # add optimization details to miscellaneous output
            if miscout is not None:
                miscout["sel_score"] = sel_score
                miscout["sel"] = sel
                miscout.update(misc) # add dict to dict

            return pgmat, sel, self.ncross, self.nprogeny

        # estimate Pareto frontier, then choose from non-dominated points.
        elif self.method == "pareto":
            # raises error
            self.pareto(
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

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
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

        G = self._calc_G(pgmat, gmat)       # get genomic relationship matrix: (n,n)

        # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
        # if we are unable to fix, then raise value error
        if not G.apply_jitter():
            raise ValueError(
                "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                "    This could be caused by lack of genetic diversity.\n"
            )

        K = G.mat_asformat("kinship")       # convert G to (1/2)G (kinship analogue): (n,n)
        C = numpy.linalg.cholesky(K).T      # cholesky decomposition of K matrix: (n,n)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,     # byte code pointer
            self.objfn_static.__globals__,  # global variables
            None,                           # new name for the function
            (C, trans, trans_kwargs),       # default values for last 3 arguments
            self.objfn_static.__closure__   # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a vectorized objective function.
        """
        # get default parameters
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        G = self._calc_G(pgmat, gmat)       # get genomic relationship matrix: (n,n)

        # to ensure we're able to perform cholesky decomposition, apply jitter if needed.
        # if we are unable to fix, then raise value error
        if not G.apply_jitter():
            raise ValueError(
                "Unable to construct objective function: Kinship matrix is not positive definite.\n"+
                "    This could be caused by lack of genetic diversity.\n"
            )

        K = G.mat_asformat("kinship")       # convert G to (1/2)G (kinship analogue): (n,n)
        C = numpy.linalg.cholesky(K).T      # cholesky decomposition of K matrix: (n,n)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (C, trans, trans_kwargs),           # default values for last 3 arguments
            self.objfn_vec_static.__closure__   # closure byte code pointer
        )

        return outfn

    def pareto(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
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
        raise RuntimeError("BinaryMinimumMeanGenomicRelationshipSelection is single objective")

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel: numpy.ndarray, C: numpy.ndarray, trans: Callable, kwargs: dict):
        """
        Score a parent contribution vector according to its mean genomic relationship.

        Parameters
        ----------
        sel : numpy.ndarray
            A parent contribution vector of shape ``(n,)`` and floating dtype.

            Where:

            - ``n`` is the number of individuals.
        C : numpy.ndarray
            An upper triangle matrix of shape (n,n) resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

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
        mgr : numpy.ndarray
            A matrix of shape (1,) if ``trans`` is ``None``.

            The first index in the array is the mean genomic relationship:

            .. math::
                MGR = || \\textbf{C} \\textbf{(sel)} ||_2

            Other indices are the mean expected trait values for the other ``t``
            traits. Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate MGR
        # (n,n)[:,(k,)] -> (n,k)
        # (1/k) * (n,k).sum(1) -> (n,)
        Cx = (1.0 / sel.shape[0]) * C[:,sel].sum(1)

        # norm2( (n,) ) -> scalar
        mgr = numpy.linalg.norm(Cx, ord = 2)

        # [scalar] -> (1,)
        mgr = numpy.array([mgr])
        
        # apply transformations if needed
        if trans:
            mgr = trans(mgr, **kwargs)
        
        return mgr

    @staticmethod
    def objfn_vec_static(sel: numpy.ndarray, C: numpy.ndarray, trans: Callable, kwargs: dict):
        """
        Score a parent contribution vector according to its mean genomic relationship.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of configurations to score.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            ``sel`` cannot be ``None``.
        C : numpy.ndarray
            An upper triangle matrix of shape (n,n) resulting from a Cholesky 
            decomposition of half a genomic relationship matrix: (1/2)G = C'C.

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
        mgr : numpy.ndarray
            A matrix of shape (j,1) if ``trans`` is ``None``.

            The first index in the array is the mean genomic relationship:

            .. math::
                MGR = 1 - || \\textbf{C} \\textbf{(sel)} ||_2

            Other indices are the mean expected trait values for the other ``t``
            traits. Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate MEH
        # (n,n)[:,(j,k)] -> (n,j,k)
        # (n,j,k).sum(2) -> (n,j)
        Cx = (1.0 / sel.shape[1]) * C[:,sel].sum(2)
        
        # norm2( (n,j), axis=0 ) -> (j,)
        # (j,)[:,None] -> (j,1)
        mgr = numpy.linalg.norm(Cx, ord = 2, axis = 0)[:,None]

        # apply transformations if needed
        if trans:
            mgr = trans(mgr, **kwargs)
        
        return mgr
