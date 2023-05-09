"""
Module implementing selection protocols for conventional phenotypic selection.
"""

from typing import Callable, Optional, Union
import numpy
import types
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm, check_is_OptimizationAlgorithm
from pybrops.opt.algo.NSGA2SetGeneticAlgorithm import NSGA2SetGeneticAlgorithm
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_gt
from pybrops.core.error import check_is_str
from pybrops.core.error import check_is_dict
from pybrops.core.error import check_is_callable
from pybrops.core.error import check_is_Generator_or_RandomState
from pybrops.core.random import global_prng

class ConventionalPhenotypicSelection(UnconstrainedSelectionProtocol):
    """
    Class implementing selection protocols for conventional phenotypic selection.

    # TODO: add formulae for methodology.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            nparent: int, 
            ncross: int, 
            nprogeny: int,
            method: str = "single",
            objfn_trans = None, 
            objfn_trans_kwargs = None, 
            objfn_wt = 1.0,
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = 1.0,
            rng = global_prng, 
            moalgo: Optional[OptimizationAlgorithm] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for Conventional Phenotypic Selection (CPS)

        Parameters
        ----------
        nparent : int
        ncross : int
        nprogeny : int
        objfn_trans : function, callable, None
        objfn_trans_kwargs : dict, None
        objfn_wt : float, numpy.ndarray
        ndset_trans : function, callable, None
        ndset_trans_kwargs : dict, None
        ndset_wt : float
        rng : numpy.random.Generator, numpy.random.RandomState
        """
        super(ConventionalPhenotypicSelection, self).__init__(**kwargs)

        # variable assignment
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = rng
        self.moalgo = moalgo

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nparent(self) -> int:
        """Number of parents to select."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set number of parents to select."""
        check_is_int(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value
    @nparent.deleter
    def nparent(self) -> None:
        """Delete number of parents to select."""
        del self._nparent

    @property
    def ncross(self) -> int:
        """Number of crosses per configuration."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
        """Set number of crosses per configuration."""
        check_is_int(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value
    @ncross.deleter
    def ncross(self) -> None:
        """Delete number of crosses per configuration."""
        del self._ncross

    @property
    def nprogeny(self) -> int:
        """Number of progeny to derive from each cross configuration."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set number of progeny to derive from each cross configuration."""
        check_is_int(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value
    @nprogeny.deleter
    def nprogeny(self) -> None:
        """Delete number of progeny to derive from each cross configuration."""
        del self._nprogeny

    @property
    def method(self) -> str:
        """Selection method."""
        return self._method
    @method.setter
    def method(self, value: str) -> None:
        """Set selection method."""
        check_is_str(value, "method")       # must be string
        value = value.lower()               # convert to lowercase
        options = ("single", "pareto")      # method options
        # if not method supported raise ValueError
        if value not in options:
            raise ValueError("Unsupported 'method'. Options are: " + ", ".join(map(str, options)))
        self._method = value
    @method.deleter
    def method(self) -> None:
        """Delete selection method."""
        del self._method

    @property
    def objfn_trans(self) -> Union[Callable,None]:
        """Objective function transformation function."""
        return self._objfn_trans
    @objfn_trans.setter
    def objfn_trans(self, value: Union[Callable,None]) -> None:
        """Set objective function transformation function."""
        if value is not None:                       # if given object
            check_is_callable(value, "objfn_trans") # must be callable
        self._objfn_trans = value
    @objfn_trans.deleter
    def objfn_trans(self) -> None:
        """Delete objective function transformation function."""
        del self._objfn_trans

    @property
    def objfn_trans_kwargs(self) -> dict:
        """Objective function transformation function keyword arguments."""
        return self._objfn_trans_kwargs
    @objfn_trans_kwargs.setter
    def objfn_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set objective function transformation function keyword arguments."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "objfn_trans_kwargs")  # check is dict
        self._objfn_trans_kwargs = value
    @objfn_trans_kwargs.deleter
    def objfn_trans_kwargs(self) -> None:
        """Delete objective function transformation function keyword arguments."""
        del self._objfn_trans_kwargs

    @property
    def objfn_wt(self) -> Union[float,numpy.ndarray]:
        """Objective function weights."""
        return self._objfn_wt
    @objfn_wt.setter
    def objfn_wt(self, value: Union[float,numpy.ndarray]) -> None:
        """Set objective function weights."""
        self._objfn_wt = value
    @objfn_wt.deleter
    def objfn_wt(self) -> None:
        """Delete objective function weights."""
        del self._objfn_wt

    @property
    def ndset_trans(self) -> Union[Callable,None]:
        """Nondominated set transformation function."""
        return self._ndset_trans
    @ndset_trans.setter
    def ndset_trans(self, value: Union[Callable,None]) -> None:
        """Set nondominated set transformation function."""
        if value is not None:                       # if given object
            check_is_callable(value, "ndset_trans") # must be callable
        self._ndset_trans = value
    @ndset_trans.deleter
    def ndset_trans(self) -> None:
        """Delete nondominated set transformation function."""
        del self._ndset_trans

    @property
    def ndset_trans_kwargs(self) -> dict:
        """Nondominated set transformation function keyword arguments."""
        return self._ndset_trans_kwargs
    @ndset_trans_kwargs.setter
    def ndset_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set nondominated set transformation function keyword arguments."""
        if value is None:                           # if given None
            value = {}                              # set default to empty dict
        check_is_dict(value, "ndset_trans_kwargs")  # check is dict
        self._ndset_trans_kwargs = value
    @ndset_trans_kwargs.deleter
    def ndset_trans_kwargs(self) -> None:
        """Delete nondominated set transformation function keyword arguments."""
        del self._ndset_trans_kwargs

    @property
    def ndset_wt(self) -> Union[float,numpy.ndarray]:
        """Nondominated set weights."""
        return self._ndset_wt
    @ndset_wt.setter
    def ndset_wt(self, value: Union[float,numpy.ndarray]) -> None:
        """Set nondominated set weights."""
        self._ndset_wt = value
    @ndset_wt.deleter
    def ndset_wt(self) -> None:
        """Delete nondominated set weights."""
        del self._ndset_wt

    @property
    def rng(self) -> Union[numpy.random.Generator,numpy.random.RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[numpy.random.Generator,numpy.random.RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator source."""
        del self._rng

    @property
    def moalgo(self) -> OptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[OptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            value = NSGA2SetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                mu = 100,       # number of parents in population
                lamb = 100,     # number of progeny to produce
                M = 1.5,        # algorithm crossover genetic map length
                rng = self.rng  # PRNG source
            )
        check_is_OptimizationAlgorithm(value, "moalgo")
        self._moalgo = value
    @moalgo.deleter
    def moalgo(self) -> None:
        """Delete multi-objective opimization algorithm."""
        del self._moalgo

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
            Genotypes (unphased most likely)
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
        method : str
            Options: "single", "pareto"
        nparent : int
        ncross : int
        nprogeny : int
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
        # single objective method: objfn_trans returns a single value for each
        # selection configuration
        if self.method == "single":
            # get vectorized objective function
            objfn = self.objfn(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max,
                **kwargs
            )

            # get all EBVs for each individual
            # (n,)
            ebv = [objfn([i]) for i in range(bvmat.ntaxa)]

            # convert to numpy.ndarray
            ebv = numpy.array(ebv)

            # multiply the objectives by objfn_wt to transform to maximizing function
            # (n,) * scalar -> (n,)
            ebv = ebv * self.objfn_wt

            # get indices of top nparent EBVs
            sel = ebv.argsort()[::-1][:self.nparent]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

            # get GEBVs for reference
            if miscout is not None:
                miscout["ebv"] = ebv

            return pgmat, sel, self.ncross, self.nprogeny

        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
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

        # raise error since method is not supported
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")

    def objfn(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return an objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Not used by this function.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Used by this function. Input breeding value matrix.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A selection objective function for the specified problem.
        """
        # get default parameters if any are None
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # get pointers to raw numpy.ndarray matrices
        mat = bvmat.mat     # (n,t) get breeding value matrix

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,     # byte code pointer
            self.objfn_static.__globals__,  # global variables
            None,                           # new name for the function
            (mat, trans, trans_kwargs),     # default values for arguments
            self.objfn_static.__closure__   # closure byte code pointer
        )

        return outfn

    def objfn_vec(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, **kwargs: dict):
        """
        Return a vectorized objective function for the provided datasets.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Not used by this function.
        gmat : GenotypeMatrix
            Not used by this function.
        ptdf : PhenotypeDataFrame
            Not used by this function.
        bvmat : BreedingValueMatrix
            Used by this function. Input breeding value matrix.
        gpmod : LinearGenomicModel
            Linear genomic prediction model.

        Returns
        -------
        outfn : function
            A vectorized selection objective function for the specified problem.
        """
        # get default parameters if any are None
        trans = self.objfn_trans
        trans_kwargs = self.objfn_trans_kwargs

        # get pointers to raw numpy.ndarray matrices
        mat = bvmat.mat     # (n,t) get breeding value matrix

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (mat, trans, trans_kwargs),         # default values for arguments
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
        miscout : dict, None, default = None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple containing two objects ``(frontier, sel_config)``.

            Where:

            - frontier is a ``numpy.ndarray`` of shape ``(q,v)`` containing
              Pareto frontier points.
            - sel_config is a ``numpy.ndarray`` of shape ``(q,k)`` containing
              parent selection decisions for each corresponding point in the
              Pareto frontier.

            Where:

            - ``q`` is the number of points in the frontier.
            - ``v`` is the number of objectives for the frontier.
            - ``k`` is the number of search space decision variables.
        """
        # get number of taxa
        ntaxa = bvmat.ntaxa

        # create objective function
        objfn = self.objfn(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max,
            **kwargs
        )

        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,                  # objective function
            k = self.nparent,               # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = self.objfn_wt        # weights to apply to each objective
        )

        # handle miscellaneous output
        if miscout is not None:     # if miscout is provided
            miscout.update(misc)    # add 'misc' to 'miscout', overwriting as needed

        return frontier, sel_config

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    @staticmethod
    def objfn_static(sel, mat, trans, kwargs):
        """
        Score a selection configuration based on its breeding values
        (Conventional Phenotype Selection; CPS).

        CPS selects the ``q`` individuals with the largest EBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is ``None``, use all individuals.
        mat : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cps : numpy.ndarray
            A EBV matrix of shape ``(t,)`` if ``trans`` is ``None``.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``t`` is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            sel = slice(None)

        # CPS calculation explanation
        # Step 1: (n,t) -> (k,t)        # select individuals
        # Step 2: (k,t).sum(0) -> (t,)  # sum across all individuals
        cps = mat[sel,:].sum(0)

        # apply transformations
        # (t,) ---trans---> (?,)
        if trans:
            cps = trans(cps, **kwargs)

        return cps

    @staticmethod
    def objfn_vec_static(sel, mat, trans, kwargs):
        """
        Score a selection configuration based on its breeding values
        (Conventional Phenotype Selection; CPS).

        CPS selects the ``q`` individuals with the largest EBVs.

        Parameters
        ----------
        sel : numpy.ndarray, None
            A selection indices matrix of shape ``(j,k)``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
            If ``sel`` is ``None``, score each individual separately: ``(n,1)``
        mat : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        trans : function or callable
            A transformation operator to alter the output.
            Function must adhere to the following standard:

            - Must accept a single ``numpy.ndarray`` argument.
            - Must return a single object, whether scalar or ``numpy.ndarray``.
        kwargs : dict
            Dictionary of keyword arguments to pass to ``trans`` function.

        Returns
        -------
        cps : numpy.ndarray
            A GEBV matrix of shape ``(j,t)`` if ``trans`` is None.
            Otherwise, of shape specified by ``trans``.

            Where:

            - ``j`` is the number of selection configurations.
            - ``t`` is the number of traits.
        """
        # if sel is None, slice all individuals
        if sel is None:
            n = mat.shape[0]
            sel = numpy.arange(n).reshape(n,1)

        # CPS calculation explanation
        # Step 1: (n,t)[(j,k),:] -> (j,k,t) # select individuals
        # Step 2: (j,k,t).sum(1) -> (j,t)   # sum across all individuals
        cps = mat[sel,:].sum(1)

        # apply transformations
        # (j,t) ---trans---> (?,?)
        if trans:
            cps = trans(cps, **kwargs)

        return cps
