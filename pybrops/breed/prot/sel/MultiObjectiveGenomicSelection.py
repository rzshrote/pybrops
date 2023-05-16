"""
Module implementing selection protocols for Optimal Population Value selection.
"""

__all__ = [
    "MultiObjectiveGenomicBaseSelection",
    "MultiObjectiveGenomicSubsetSelection"
]

from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Tuple, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.MultiObjectiveGenomicSelectionProblem import MultiObjectiveGenomicSubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel, check_is_AdditiveLinearGenomicModel
from pybrops.opt.algo.ConstrainedNSGA2SubsetGeneticAlgorithm import ConstrainedNSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm, check_is_ConstrainedOptimizationAlgorithm
from pybrops.core.random.prng import global_prng
from pybrops.opt.algo.ConstrainedSteepestDescentSubsetHillClimber import ConstrainedSteepestDescentSubsetHillClimber
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class MultiObjectiveGenomicSelection(SelectionProtocol,metaclass=ABCMeta):
    """
    Semi-abstract class for Optimal Population Value (OPV) Selection with constraints.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
            weight: Union[numpy.ndarray,Callable],
            target: Union[numpy.ndarray,Callable],
            method: str,
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[ConstrainedOptimizationAlgorithm] = None,
            moalgo: Optional[ConstrainedOptimizationAlgorithm] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for MultiObjectiveGenomicSelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(MultiObjectiveGenomicSelection, self).__init__(
            method = method,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            ndset_wt = ndset_wt,
            ndset_trans = ndset_trans,
            ndset_trans_kwargs = ndset_trans_kwargs
            **kwargs
        )
        # order dependent assignments
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.weight = weight
        self.target = target
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################ Object Properties #############################
    @property
    def nparent(self) -> Integral:
        """Number of parents to select."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set number of parents to select."""
        check_is_Integral(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value

    @property
    def ncross(self) -> Integral:
        """Number of crosses per configuration."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: Integral) -> None:
        """Set number of crosses per configuration."""
        check_is_Integral(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value

    @property
    def nprogeny(self) -> Integral:
        """Number of progeny to derive from each cross configuration."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: Integral) -> None:
        """Set number of progeny to derive from each cross configuration."""
        check_is_Integral(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value

    @property
    def weight(self) -> Union[numpy.ndarray,Callable]:
        """Allele weights."""
        return self._weight
    @weight.setter
    def weight(self, value: Union[numpy.ndarray,Callable]) -> None:
        """Set allele weights."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "weight", 2)
        elif callable(value):
            pass
        else:
            raise TypeError("variable 'weight' must be a callable function or numpy.ndarray")
        self._weight = value
    
    @property
    def target(self) -> Union[numpy.ndarray,Callable]:
        """Target allele frequency."""
        return self._target
    @target.setter
    def target(self, value: Union[numpy.ndarray,Callable]) -> None:
        """Set target allele frequency."""
        if isinstance(value, numpy.ndarray):
            check_ndarray_ndim(value, "target", 2)
        elif callable(value):
            pass
        else:
            raise TypeError("variable 'target' must be a callable function or numpy.ndarray")
        self._target = value

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value

    @property
    @abstractmethod
    def soalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        raise NotImplementedError("method is abstract")
    @soalgo.setter
    @abstractmethod
    def soalgo(self, value: ConstrainedOptimizationAlgorithm) -> None:
        """Set single-objective optimization algorithm."""
        raise NotImplementedError("method is abstract")

    @property
    @abstractmethod
    def moalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        raise NotImplementedError("method is abstract")
    @moalgo.setter
    @abstractmethod
    def moalgo(self, value: ConstrainedOptimizationAlgorithm) -> None:
        """Set multi-objective opimization algorithm."""
        raise NotImplementedError("method is abstract")

    ######################### Private Object Methods ###########################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    # leave problem() abstract

    ############## Pareto Frontier Functions ###############
    def pareto(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: AdditiveLinearGenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: dict = None, 
            **kwargs: dict
        ) -> Tuple[numpy.ndarray,numpy.ndarray]:
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
        # construct the optimization problem
        prob = self.problem(
            pgmat = pgmat,
            gmat = gmat,
            ptdf = ptdf,
            bvmat = bvmat,
            gpmod = gpmod,
            t_cur = t_cur,
            t_max = t_max
        )

        # optimize the problem
        soln = self.moalgo.minimize(
            prob = prob,
            miscout = miscout
        )

        if miscout is not None:
            miscout["soln"] = soln

        frontier = soln.soln_obj
        sel_config = soln.soln_decn

        return frontier, sel_config

    ################# Selection Functions ##################
    def select(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: AdditiveLinearGenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: dict = None, 
            **kwargs: dict
        ) -> Tuple[PhasedGenotypeMatrix,numpy.ndarray,numpy.ndarray,numpy.ndarray]:
        """
        Select parents individuals for breeding.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes
        gmat : GenotypeMatrix
            Genotype matrix from which to calculate genomic relationships.
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
        # check inputs
        check_is_PhasedGenotypeMatrix(pgmat, "pgmat")
        check_is_AdditiveLinearGenomicModel(gpmod, "gpmod")

        # Solve problem using a single objective method
        if self.method == "single":
            # create optimization problem
            prob = self.problem(
                pgmat = pgmat,
                gmat = gmat,
                ptdf = ptdf,
                bvmat = bvmat,
                gpmod = gpmod,
                t_cur = t_cur,
                t_max = t_max
            )

            # miscellaneous output
            misc = {}

            # use single-objective optimization to get a solution.
            soln = self.soalgo.minimize(
                prob = prob,
                miscout = miscout
            )

            # get first and only solution
            obj = soln.soln_obj[0]
            sel = soln.soln_decn[0]

            # add optimization details to miscellaneous output
            if miscout is not None:
                miscout["obj"] = obj
                miscout["sel"] = sel
                miscout.update(misc) # add dict to dict

            return pgmat, sel, self.ncross, self.nprogeny

        # estimate Pareto frontier, then choose from non-dominated points.
        elif self.method == "pareto":
            # raises error
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

class MultiObjectiveGenomicSubsetSelection(MultiObjectiveGenomicSelection):
    """
    Class defining Optimal Haploid Value (OHV) Selection for subset search spaces.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################
    @property
    def soalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[ConstrainedOptimizationAlgorithm,None]) -> None:
        """Set single-objective optimization algorithm."""
        if value is None:
            # construct default hillclimber
            value = ConstrainedSteepestDescentSubsetHillClimber(self.rng)
        check_is_ConstrainedOptimizationAlgorithm(value, "soalgo")
        self._soalgo = value

    @property
    def moalgo(self) -> ConstrainedOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[ConstrainedOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = ConstrainedNSGA2SubsetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100, # number of parents in population
                rng = self.rng  # PRNG source
            )
        check_is_ConstrainedOptimizationAlgorithm(value, "moalgo")
        self._moalgo = value

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: AdditiveLinearGenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            **kwargs: dict
        ) -> SelectionProblem:
        """
        Create an optimization problem definition using provided inputs.

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
        out : SelectionProblem
            An optimization problem definition.
        """
        # get decision space parameters
        ntaxa = gmat.ntaxa
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = MultiObjectiveGenomicSubsetSelectionProblem.from_object(
            gmat = gmat,
            weight = self.weight,
            target = self.target,
            gpmod = gpmod,
            ndecn = self.nparent,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = self.nobj,
            obj_wt = self.obj_wt,
            obj_trans = self.obj_trans,
            obj_trans_kwargs = self.obj_trans_kwargs,
            nineqcv = self.nineqcv,
            ineqcv_wt = self.ineqcv_wt,
            ineqcv_trans = self.ineqcv_trans,
            ineqcv_trans_kwargs = self.ineqcv_trans_kwargs,
            neqcv = self.neqcv,
            eqcv_wt = self.eqcv_wt,
            eqcv_trans = self.eqcv_trans,
            eqcv_trans_kwargs = self.eqcv_trans_kwargs
        )

        return prob

    ############## Pareto Frontier Functions ###############
    # inherit pareto() implementation

    ################# Selection Functions ##################
    # inherit select() implementation
