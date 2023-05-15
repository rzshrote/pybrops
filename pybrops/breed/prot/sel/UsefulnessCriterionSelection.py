"""
Module implementing selection protocols for Usefulness Criterion (UC) selection.
"""

__all__ = [
    "UsefulnessCriterionBaseSelection",
    "UsefulnessCriterionBinarySelection",
    "UsefulnessCriterionIntegerSelection",
    "UsefulnessCriterionRealSelection",
    "UsefulnessCriterionSubsetSelection"
]

from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
import numpy
import scipy.stats
from numpy.random import Generator, RandomState
from typing import Optional, Union
from typing import Callable
from pybrops.breed.prot.sel.prob.UsefulnessCriterionSelectionProblem import UsefulnessCriterionBinarySelectionProblem, UsefulnessCriterionIntegerSelectionProblem, UsefulnessCriterionRealSelectionProblem, UsefulnessCriterionSubsetSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.model.vmat.fcty.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory, check_is_GeneticVarianceMatrixFactory
from pybrops.opt.algo.ConstrainedNSGA2SubsetGeneticAlgorithm import ConstrainedNSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.ConstrainedSteepestDescentSubsetHillClimber import ConstrainedSteepestDescentSubsetHillClimber
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm, check_is_ConstrainedOptimizationAlgorithm
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Integral, check_is_Real, check_is_bool
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq, check_is_in_interval, check_is_lt
from pybrops.core.random.prng import global_prng
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix, check_is_BreedingValueMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix, check_is_GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix, check_is_PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class UsefulnessCriterionBaseSelection(SelectionProtocol,metaclass=ABCMeta):
    """
    Semi-abstract class for Optimal Haploid Value (OHV) Selection with constraints.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            nconfig: Integral, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral,
            upper_percentile: Real,
            vmatfcty: GeneticVarianceMatrixFactory,
            gmapfn: GeneticMapFunction,
            unique_parents: bool, 
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
        Constructor for ConventionalGenomicSelection.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(UsefulnessCriterionBaseSelection, self).__init__(
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
        self.nconfig = nconfig
        self.nparent = nparent
        self.ncross = ncross
        self.nprogeny = nprogeny
        self.nself = nself
        self.upper_percentile = upper_percentile
        self.vmatfcty = vmatfcty
        self.gmapfn = gmapfn
        self.unique_parents = unique_parents
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################ Object Properties #############################
    @property
    def nconfig(self) -> Integral:
        """Number of cross configurations to consider."""
        return self._nconfig
    @nconfig.setter
    def nconfig(self, value: Integral) -> None:
        """Set number of cross configurations to consider."""
        check_is_Integral(value, "nconfig")      # must be int
        check_is_gt(value, "nconfig", 0)    # int must be >0
        self._nconfig = value

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
    def nself(self) -> Integral:
        """Description for property nself."""
        return self._nself
    @nself.setter
    def nself(self, value: Integral) -> None:
        """Set data for property nself."""
        check_is_Integral(value, "nself")     # must be int
        check_is_gteq(value, "nself", 0)   # int must be >=0
        self._nself = value

    @property
    def upper_percentile(self) -> Real:
        """Description for property upper_percentile."""
        return self._upper_percentile
    @upper_percentile.setter
    def upper_percentile(self, value: Real) -> None:
        """Set data for property upper_percentile."""
        check_is_Real(value, "upper_percentile")  # must be a number
        check_is_gt(value, "upper_percentile", 0.0)   # number must be >0
        check_is_lt(value, "upper_percentile", 1.0)   # number must be <1
        self._upper_percentile = value
        # set the selection intensity
        self._selection_intensity = scipy.stats.norm.pdf(scipy.stats.norm.ppf(1.0 - self._upper_percentile)) / self._upper_percentile
    
    @property
    def selection_intensity(self) -> Real:
        """Description for property selection_intensity."""
        return self._selection_intensity
    @selection_intensity.setter
    def selection_intensity(self, value: Real) -> None:
        """Set data for property selection_intensity."""
        error_readonly("selection_intensity")

    @property
    def vmatfcty(self) -> GeneticVarianceMatrixFactory:
        """Description for property vmatfcty."""
        return self._vmatfcty
    @vmatfcty.setter
    def vmatfcty(self, value: GeneticVarianceMatrixFactory) -> None:
        """Set data for property vmatfcty."""
        check_is_GeneticVarianceMatrixFactory(value, "vmatfcty")
        self._vmatfcty = value
    
    @property
    def gmapfn(self) -> GeneticMapFunction:
        """Description for property gmapfn."""
        return self._gmapfn
    @gmapfn.setter
    def gmapfn(self, value: GeneticMapFunction) -> None:
        """Set data for property gmapfn."""
        self._gmapfn = value

    @property
    def unique_parents(self) -> bool:
        """Whether parents should be unique."""
        return self._unique_parents
    @unique_parents.setter
    def unique_parents(self, value: bool) -> None:
        """Set whether parents should be unique."""
        check_is_bool(value, "unique_parents")
        self._unique_parents = value

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
    def select(self, pgmat, gmat, ptdf, bvmat, gpmod, t_cur, t_max, miscout = None, **kwargs: dict):
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
        check_is_GenotypeMatrix(gmat, "gmat")
        check_is_BreedingValueMatrix(bvmat, "bvmat")
        check_is_Real(t_cur, "t_cur")
        check_is_Real(t_max, "t_max")

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

class UsefulnessCriterionSubsetSelection(UsefulnessCriterionBaseSelection):
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
            gpmod: GenomicModel, 
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
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionSubsetSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space = numpy.arange(len(xmap))
        decn_space_lower = numpy.repeat(0, self.nconfig)
        decn_space_upper = numpy.repeat(len(xmap)-1, self.nconfig)

        # construct problem
        prob = UsefulnessCriterionSubsetSelectionProblem.from_object(
            nparent = self.nparent, 
            ncross = self.ncross, 
            nprogeny = self.nprogeny, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = self.nconfig,
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

class UsefulnessCriterionRealSelection(UsefulnessCriterionBaseSelection):
    """
    Class defining Optimal Haploid Value (OHV) Selection for real search spaces.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionRealSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, len(xmap))
        decn_space_upper = numpy.repeat(1.0, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = UsefulnessCriterionRealSelectionProblem.from_object(
            nparent = self.nparent, 
            ncross = self.ncross, 
            nprogeny = self.nprogeny, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = len(xmap),
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

class UsefulnessCriterionIntegerSelection(UsefulnessCriterionBaseSelection):
    """
    Class defining Optimal Haploid Value (OHV) Selection for a integer search spaces.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionIntegerSelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(self.nconfig * self.nparent * self.ncross, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = UsefulnessCriterionIntegerSelectionProblem.from_object(
            nparent = self.nparent, 
            ncross = self.ncross, 
            nprogeny = self.nprogeny, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = len(xmap),
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

class UsefulnessCriterionBinarySelection(UsefulnessCriterionBaseSelection):
    """
    Class defining Optimal Haploid Value (OHV) Selection for a binary search spaces.
    """

    ########################## Special Object Methods ##########################
    # inherit __init__() implementation

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    def problem(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
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
        # get the cross map (inefficient)
        xmap = UsefulnessCriterionBinarySelectionProblem._calc_xmap(
            pgmat.ntaxa,
            self.nparent,
            self.unique_parents
        )

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, len(xmap))
        decn_space_upper = numpy.repeat(1, len(xmap))
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = UsefulnessCriterionBinarySelectionProblem.from_object(
            nparent = self.nparent, 
            ncross = self.ncross, 
            nprogeny = self.nprogeny, 
            nself = self.nself,
            upper_percentile = self.upper_percentile,
            vmatfcty = self.vmatfcty,
            gmapfn = self.gmapfn,
            unique_parents = self.unique_parents, 
            pgmat = pgmat,
            gpmod = gpmod,
            ndecn = len(xmap),
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
