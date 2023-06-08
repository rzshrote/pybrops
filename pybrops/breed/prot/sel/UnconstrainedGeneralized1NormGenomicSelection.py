"""
Module implementing selection protocols for general 1-norm genomic selection
"""

import types
import numpy
from typing import Union
from typing import Callable

from pybrops.opt.algo.UnconstrainedNSGA2SetGeneticAlgorithm import UnconstrainedNSGA2SetGeneticAlgorithm
from pybrops.opt.algo.UnconstrainedOptimizationAlgorithm import UnconstrainedOptimizationAlgorithm, check_is_OptimizationAlgorithm
from pybrops.opt.algo.UnconstrainedSteepestAscentSetHillClimber import UnconstrainedSteepestAscentSetHillClimber
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.breed.prot.sel.targetfn import target_positive
from pybrops.breed.prot.sel.weightfn import weight_absolute
from pybrops.core.error.error_attr_python import check_is_callable
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_type_python import check_is_int
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_python import check_isinstance
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class Generalized1NormGenomicSelection(UnconstrainedSelectionProtocol):
    """
    docstring for Generalized1NormGenomicSelection.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            nparent: int,
            ncross: int,
            nprogeny: int,
            weight: Union[numpy.ndarray,Callable] = weight_absolute,
            target: Union[numpy.ndarray,Callable] = target_positive,
            method: str = "single",
            objfn_trans = None,
            objfn_trans_kwargs = None, 
            objfn_wt = -1.0,
            ndset_trans = None, 
            ndset_trans_kwargs = None, 
            ndset_wt = -1.0,
            rng = global_prng, 
            soalgo = None, 
            moalgo = None,
            **kwargs
        ):
        """
        Constructor for Generalized1NormGenomicSelection.
        
        Parameters
        ----------
        nparent : int
            Number of parents to select.
        """
        super(Generalized1NormGenomicSelection, self).__init__(**kwargs)

        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.weight = weight
        self.target = target
        self.method = method
        self.objfn_trans = objfn_trans
        self.objfn_trans_kwargs = objfn_trans_kwargs
        self.objfn_wt = objfn_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.ndset_wt = ndset_wt
        self.rng = rng
        self.soalgo = soalgo
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
    def weight(self) -> Union[numpy.ndarray,Callable]:
        """Marker weights or marker weight function."""
        return self._weight
    @weight.setter
    def weight(self, value: Union[numpy.ndarray,Callable]) -> None:
        """Set marker weights or marker weight function."""
        check_isinstance(value, "weight", (numpy.ndarray, Callable))
        self._weight = value
    @weight.deleter
    def weight(self) -> None:
        """Delete marker weights or marker weight function."""
        del self._weight

    @property
    def target(self) -> Union[numpy.ndarray,Callable]:
        """Allele frequency targets or allele frequency target function."""
        return self._target
    @target.setter
    def target(self, value: Union[numpy.ndarray,Callable]) -> None:
        """Set allele frequency targets or allele frequency target function."""
        check_isinstance(value, "target", (numpy.ndarray, Callable))
        self._target = value
    @target.deleter
    def target(self) -> None:
        """Delete allele frequency targets or allele frequency target function."""
        del self._target

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
    def soalgo(self) -> UnconstrainedOptimizationAlgorithm:
        """Single objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[UnconstrainedOptimizationAlgorithm,None]) -> None:
        """Set single objective optimization algorithm."""
        if value is None:
            value = UnconstrainedSteepestAscentSetHillClimber(rng = self.rng)
        check_is_OptimizationAlgorithm(value, "soalgo")
        self._soalgo = value
    @soalgo.deleter
    def soalgo(self) -> None:
        """Delete single objective optimization algorithm."""
        del self._soalgo

    @property
    def moalgo(self) -> UnconstrainedOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[UnconstrainedOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            value = UnconstrainedNSGA2SetGeneticAlgorithm(
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
    ########################## Private Object Methods ##########################
    ############################################################################
    def _calc_mkrwt(self, gpmod: AdditiveLinearGenomicModel):
        if callable(self.weight):
            return self.weight(gpmod.u_a)
        elif isinstance(self.weight, numpy.ndarray):
            return self.weight
        else:
            raise TypeError("variable 'weight' must be a callable function or numpy.ndarray")
    
    def _calc_tfreq(self, gpmod: AdditiveLinearGenomicModel):
        if callable(self.target):
            return self.target(gpmod.u_a)
        elif isinstance(self.target, numpy.ndarray):
            return self.target
        else:
            raise TypeError("variable 'target' must be a callable function or numpy.ndarray")

    def _calc_V(self, gmat: GenotypeMatrix, gpmod: AdditiveLinearGenomicModel):
        # get marker weights and target frequencies
        mkrwt = self._calc_mkrwt(gpmod) # (p,t)
        tfreq = self._calc_tfreq(gpmod) # (p,t)

        # make sure mkrwt and tfreq have the same dimensions
        if mkrwt.shape != tfreq.shape:
            raise ValueError("marker weights and target allele frequencies do not have the same shape")

        # get number of traits and taxa
        nvrnt = mkrwt.shape[0]
        ntrait = mkrwt.shape[1]
        ntaxa = gmat.ntaxa
        ploidy = float(gmat.ploidy)

        # allocate a tensor for storing distance matrices
        # (n,p,t)
        V = numpy.empty((ntaxa,nvrnt,ntrait), dtype = "float64")

        # get the genotype matrix as {0,1,2,...}
        # divide the genotypes by the number of phases
        # (n,p)
        Z = (1.0 / ploidy) * gmat.tacount()

        # calculate a distance matrix for each trait
        for trait in range(ntrait):
            # (n,p) = (1,p) * ( (n,p) - (1,p) )
            V[:,:,trait] = mkrwt[None,:,trait] * (Z - tfreq[None,:,trait])

        # return distance value tensor
        # (n,p,t)
        return V

    ############################## Object Methods ##############################
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
        # Solve problem using a single objective method
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

        V = self._calc_V(gmat, gpmod)       # get Cholesky decomposition of genomic relationship matrix: (t,n,n)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_static.__code__,     # byte code pointer
            self.objfn_static.__globals__,  # global variables
            None,                           # new name for the function
            (V, trans, trans_kwargs),       # default values for last 3 arguments
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

        V = self._calc_V(gmat, gpmod)       # get Cholesky decomposition of genomic relationship matrix: (t,n,n)

        # copy objective function and modify default values
        # this avoids using functools.partial and reduces function execution time.
        outfn = types.FunctionType(
            self.objfn_vec_static.__code__,     # byte code pointer
            self.objfn_vec_static.__globals__,  # global variables
            None,                               # new name for the function
            (V, trans, trans_kwargs),           # default values for last 3 arguments
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
        # get selection parameters
        nparent = self.nparent
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

        # use multi-objective optimization to approximate Pareto front.
        frontier, sel_config, misc = self.moalgo.optimize(
            objfn = objfn,                  # objective function
            k = nparent,                    # vector length to optimize (sspace^k)
            sspace = numpy.arange(ntaxa),   # search space options
            objfn_wt = objfn_wt,            # weights to apply to each objective
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
    def objfn_static(sel: numpy.ndarray, V: numpy.ndarray, trans: Callable, kwargs: dict):
        """
        Score a parent selection vector according to its distance from a utopian point.
        The goal is to minimize this function.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape ``(k,)``.

            Where:

            - ``k`` is the number of individuals to select.

            Each index indicates which individuals to select.
            Each index in ``sel`` represents a single individual's row.
        V : numpy.ndarray
            A matrix of shape ``(n,p,t)`` containing distance values of individuals' 
            alleles for each trait.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.
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
        dist : numpy.ndarray
            A matrix of shape ``(t,)`` if ``trans`` is ``None``.

            Distances for each target.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate vector
        # (n,p,t)[(k,),:] -> (k,p,t)
        # (k,p,t).sum(0) -> (p,t)
        # scalar * (p,t) -> (p,t)
        # | (p,t) | -> (p,t)
        # (p,t).sum(0) -> (t,)
        dist = numpy.absolute((1.0 / len(sel)) * V[sel,:,:].sum(0)).sum(0)

        # apply transformations if needed
        if trans:
            dist = trans(dist, **kwargs)
        
        return dist

    @staticmethod
    def objfn_vec_static(sel: numpy.ndarray, V: numpy.ndarray, trans: Callable, kwargs: dict):
        """
        Score a parent selection vector according to its distance from a utopian point.
        The goal is to minimize this function.

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
        V : numpy.ndarray
            A matrix of shape ``(n,p,t)`` containing distance values of individuals' 
            alleles for each trait.

            Where:

            - ``n`` is the number of individuals.
            - ``p`` is the number of markers.
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
        dist : numpy.ndarray
            A matrix of shape ``(j,t)`` if ``trans`` is ``None``.

            Distances for each target.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate vector
        # (n,p,t)[(j,k,),:] -> (j,k,p,t)
        # (j,k,p,t).sum(1) -> (j,p,t)
        # scalar * (j,p,t) -> (j,p,t)
        # | (j,p,t) | -> (j,p,t)
        # (j,p,t).sum(1) -> (j,t)
        dist = numpy.absolute((1.0 / len(sel)) * V[sel,:,:].sum(1)).sum(1)

        # apply transformations if needed
        if trans:
            dist = trans(dist, **kwargs)
        
        return dist
