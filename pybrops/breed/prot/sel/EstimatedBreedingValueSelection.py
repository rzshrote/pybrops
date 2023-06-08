"""
Module defining conventional genomic selection protocols.
"""

# list of all public objects in this module
__all__ = [
    "EstimatedBreedingValueBaseSelection",
    "EstimatedBreedingValueSubsetSelection",
    "EstimatedBreedingValueRealSelection",
    "EstimatedBreedingValueIntegerSelection",
    "EstimatedBreedingValueBinarySelection"
]

# imports
from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.prob.trans import trans_decnvec_sum_eq
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Callable, check_is_Integral, check_is_dict
from pybrops.core.error.error_value_python import check_is_gt, check_is_gteq
from pybrops.core.random.prng import global_prng
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.EstimatedBreedingValueSelectionProblem import EstimatedBreedingValueBinarySelectionProblem, EstimatedBreedingValueIntegerSelectionProblem, EstimatedBreedingValueRealSelectionProblem, EstimatedBreedingValueSubsetSelectionProblem
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.BinaryGeneticAlgorithm import BinaryGeneticAlgorithm
from pybrops.opt.algo.BinaryOptimizationAlgorithm import BinaryOptimizationAlgorithm, check_is_BinaryOptimizationAlgorithm
from pybrops.opt.algo.NSGA2BinaryGeneticAlgorithm import NSGA2BinaryGeneticAlgorithm
from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.opt.algo.SubsetGeneticAlgorithm import SubsetGeneticAlgorithm
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm, check_is_SubsetOptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix, check_is_BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class EstimatedBreedingValueBaseSelection(SelectionProtocol,metaclass=ABCMeta):
    """
    Semiabstract class for Conventional Genomic Selection (CGS) with constraints.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            nparent: Integral, 
            nmating: Integral, 
            nprogeny: Integral,
            nobj: Integral,
            obj_wt: Optional[numpy.ndarray] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[numpy.ndarray] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[numpy.ndarray] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            ndset_wt: Optional[Real] = None,
            ndset_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]] = None, 
            ndset_trans_kwargs: Optional[dict] = None, 
            rng: Optional[Union[Generator,RandomState]] = None, 
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for ConventionalGenomicSelection.

        Parameters
        ----------
        nobj: Integral
            Number of objectives.
            If None, the number of objectives is inferred from the input ``bvmat`` as ``bvmat.ntrait``.
        obj_wt: numpy.ndarray
            Objective function weights.
        obj_trans: Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs: dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv: Integral,
            Number of inequality constraints.
        ineqcv_wt: numpy.ndarray,
            Inequality constraint violation weights.
        ineqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an inequality constraint violation vector.
            The transformation function must be of the form: ``ineqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        ineqcv_trans_kwargs: Optional[dict],
            Keyword arguments for the latent space to inequality constraint violation space transformation function.
            If None, an empty dictionary is used.
        neqcv: Integral
            Number of equality constraints.
        eqcv_wt: numpy.ndarray
            Equality constraint violation weights.
        eqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an equality constraint violation vector.
            The transformation function must be of the form: ``eqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        eqcv_trans_kwargs: dict, None
            Keyword arguments for the latent space to equality constraint violation space transformation function.
            If None, an empty dictionary is used.
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        self.nparent = nparent
        self.nmating = nmating
        self.nprogeny = nprogeny
        self.rng = rng
        super(EstimatedBreedingValueBaseSelection, self).__init__(
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
            ndset_trans_kwargs = ndset_trans_kwargs,
            soalgo = soalgo,
            moalgo = moalgo,
            **kwargs
        )

    ############################ Object Properties #############################
    @property
    def nparent(self) -> int:
        """Number of parents to select."""
        return self._nparent
    @nparent.setter
    def nparent(self, value: int) -> None:
        """Set number of parents to select."""
        check_is_Integral(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value

    @property
    def nmating(self) -> int:
        """Number of crosses per configuration."""
        return self._ncross
    @nmating.setter
    def nmating(self, value: int) -> None:
        """Set number of crosses per configuration."""
        check_is_Integral(value, "ncross")       # must be int
        check_is_gt(value, "ncross", 0)     # int must be >0
        self._ncross = value

    @property
    def nprogeny(self) -> int:
        """Number of progeny to derive from each cross configuration."""
        return self._nprogeny
    @nprogeny.setter
    def nprogeny(self, value: int) -> None:
        """Set number of progeny to derive from each cross configuration."""
        check_is_Integral(value, "nprogeny")     # must be int
        check_is_gt(value, "nprogeny", 0)   # int must be >0
        self._nprogeny = value

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState,None]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng") # check is numpy.Generator
        self._rng = value

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    # problem() abstract function from SelectionProtocol

    ############## Pareto Frontier Functions ###############
    # pareto() abstract function from SelectionProtocol

    ################# Selection Functions ##################
    # select() abstract function from SelectionProtocol

class EstimatedBreedingValueSubsetSelection(EstimatedBreedingValueBaseSelection):
    """
    Conventional Genomic Selection in a subset search space.
    """
    ########################## Special Object Methods ##########################
    # __init__() concrete function from SelectionProtocol

    ############################ Object Properties #############################
    @property
    def soalgo(self) -> SubsetOptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[SubsetOptimizationAlgorithm,None]) -> None:
        """Set single-objective optimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = SubsetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100, # number of parents in population
                rng = self.rng  # PRNG source
            )
            # construct default hillclimber
            # value = SteepestDescentSubsetHillClimber(self.rng)
        check_is_SubsetOptimizationAlgorithm(value, "soalgo")
        self._soalgo = value

    @property
    def moalgo(self) -> SubsetOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[SubsetOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = NSGA2SubsetGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100, # number of parents in population
                rng = self.rng  # PRNG source
            )
        check_is_SubsetOptimizationAlgorithm(value, "moalgo")
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
        # get number of individuals
        ntaxa = bvmat.ntaxa

        # get decision space parameters
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = EstimatedBreedingValueSubsetSelectionProblem.from_bvmat(
            bvmat = bvmat,
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
    def mosolve(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            gmat: GenotypeMatrix, 
            ptdf: PhenotypeDataFrame, 
            bvmat: BreedingValueMatrix, 
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> tuple:
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
        miscout : dict, None
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
        # construct the problem
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
            gpmod: GenomicModel, 
            t_cur: Integral, 
            t_max: Integral, 
            miscout: Optional[dict] = None, 
            **kwargs: dict
        ) -> tuple:
        """
        Select individuals for breeding.

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
        miscout : dict, None
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
        # if the number of objectives is 1, then we use a single objective algorithm
        if self.nobj == 1:
            # construct the problem
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
            soln = self.soalgo.minimize(
                prob = prob,
                miscout = miscout
            )

            if miscout is not None:
                miscout["soln"] = soln

            # extract decision variables as a vector (ndecn,)
            sel = soln.soln_decn[0]

            return pgmat, sel, self.nmating, self.nprogeny
        # else, we use a multi-objective algorithm and apply a transformation on the non-dominated points to identify a 
        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif self.nobj > 1:
            # get the pareto frontier
            frontier, sel_config = self.mosolve(
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

            return pgmat, sel_config[ix], self.nmating, self.nprogeny
        else:
            raise ValueError("number of objectives must be greater than zero")

class EstimatedBreedingValueRealSelection(EstimatedBreedingValueBaseSelection):
    """
    Conventional Genomic Selection in a real search space.
    """
    ########################## Special Object Methods ##########################
    # __init__() concrete function from SelectionProtocol

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
        # get number of individuals
        ntaxa = bvmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, ntaxa)
        decn_space_upper = numpy.repeat(1.0, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = EstimatedBreedingValueRealSelectionProblem.from_bvmat(
            bvmat = bvmat,
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
    # use super pareto() function

    ################# Selection Functions ##################
    # use super select() function

class EstimatedBreedingValueIntegerSelection(EstimatedBreedingValueBaseSelection):
    """
    Conventional Genomic Selection in an integer search space.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    # use old init

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

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
        # get number of individuals
        ntaxa = bvmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(ntaxa, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = EstimatedBreedingValueIntegerSelectionProblem.from_bvmat(
            bvmat = bvmat,
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
    # use super pareto() function

    ################# Selection Functions ##################
    # use super select() function

class EstimatedBreedingValueBinarySelection(EstimatedBreedingValueBaseSelection):
    """
    Conventional Genomic Selection in a subset search space.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    # use BaseConventionalGenomicSelection __init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @EstimatedBreedingValueBaseSelection.neqcv.setter
    def neqcv(self, value) -> Integral:
        """Set number of equality constraint violations. If None, set to 1 for ``nparent`` summation constraint."""
        if value is None:
            value = 1
        check_is_Integral(value, "neqcv")
        check_is_gteq(value, "neqcv", 0)    # possible to have 0 equality constraints
        self._neqcv = value
    
    @EstimatedBreedingValueBaseSelection.eqcv_trans.setter
    def eqcv_trans(self, value: Union[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray],None]) -> None:
        """Set latent space to equality constraint violation transformation function. If None, set to the empty function."""
        if value is None:
            value = trans_decnvec_sum_eq
        check_is_Callable(value, "eqcv_trans")
        self._eqcv_trans = value 

    @EstimatedBreedingValueBaseSelection.eqcv_trans_kwargs.setter
    def eqcv_trans_kwargs(self, value: Union[dict,None]) -> None:
        """Set keyword arguments for the latent space to equality constraint violation transformation function. If None, set to empty dict."""
        if value is None:
            value = {"decnvec_sum": self.nparent}
        check_is_dict(value, "eqcv_trans_kwargs")
        self._eqcv_trans_kwargs = value

    @property
    def soalgo(self) -> BinaryOptimizationAlgorithm:
        """Single-objective optimization algorithm."""
        return self._soalgo
    @soalgo.setter
    def soalgo(self, value: Union[BinaryOptimizationAlgorithm,None]) -> None:
        """Set single-objective optimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = BinaryGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100, # number of parents in population
                rng = self.rng  # PRNG source
            )
            # construct default hillclimber
            # value = SteepestDescentSubsetHillClimber(self.rng)
        check_is_BinaryOptimizationAlgorithm(value, "soalgo")
        self._soalgo = value

    @property
    def moalgo(self) -> BinaryOptimizationAlgorithm:
        """Multi-objective opimization algorithm."""
        return self._moalgo
    @moalgo.setter
    def moalgo(self, value: Union[BinaryOptimizationAlgorithm,None]) -> None:
        """Set multi-objective opimization algorithm."""
        if value is None:
            # construct default multi-objective algorithm
            value = NSGA2BinaryGeneticAlgorithm(
                ngen = 250,     # number of generations to evolve
                pop_size = 100, # number of parents in population
                rng = self.rng  # PRNG source
            )
        check_is_BinaryOptimizationAlgorithm(value, "moalgo")
        self._moalgo = value

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

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
        # type checks
        check_is_BreedingValueMatrix(bvmat, "bvmat")

        # get number of individuals
        ntaxa = bvmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(1, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = EstimatedBreedingValueBinarySelectionProblem.from_bvmat(
            bvmat = bvmat,
            ndecn = ntaxa,
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
    # use BaseConventionalGenomicSelection pareto() function

    ################# Selection Functions ##################
    # use BaseConventionalGenomicSelection select() function
