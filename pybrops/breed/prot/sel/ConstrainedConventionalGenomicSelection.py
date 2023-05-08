"""
Module defining conventional genomic selection protocols.
"""

# list of all public objects in this module
__all__ = [
    "BaseConventionalGenomicSelection",
    "SubsetConventionalGenomicSelection",
    "RealConventionalGenomicSelection",
    "IntegerConventionalGenomicSelection"
]

# imports
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.ConstrainedSelectionProtocol import ConstrainedSelectionProtocol
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.random import global_prng
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.GenomicEstimatedBreedingValueSelectionProblem import GenomicEstimatedBreedingValueIntegerSelectionProblem, GenomicEstimatedBreedingValueRealSelectionProblem, GenomicEstimatedBreedingValueSubsetSelectionProblem
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.ConstrainedNSGA2SubsetGeneticAlgorithm import ConstrainedNSGA2SubsetGeneticAlgorithm
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm, check_is_ConstrainedOptimizationAlgorithm
from pybrops.opt.algo.ConstrainedSteepestDescentSubsetHillClimber import ConstrainedSteepestDescentSubsetHillClimber
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class BaseConventionalGenomicSelection(ConstrainedSelectionProtocol):
    """
    Semiabstract class for Conventional Genomic Selection (CGS) with constraints.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            nparent: Integral, 
            ncross: Integral, 
            nprogeny: Integral,
            method: str,
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
        super(BaseConventionalGenomicSelection, self).__init__(
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
        check_is_Integral(value, "nparent")      # must be int
        check_is_gt(value, "nparent", 0)    # int must be >0
        self._nparent = value

    @property
    def ncross(self) -> int:
        """Number of crosses per configuration."""
        return self._ncross
    @ncross.setter
    def ncross(self, value: int) -> None:
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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ########## Optimization Problem Construction ###########
    # leave abstract since it depends on the problem type

    ############## Pareto Frontier Functions ###############
    # implement since this is problem type agnostic
    def pareto(
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
    # implement since this is problem type agnostic
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
        # single-objective method: objfn_trans returns a single value for each
        # selection configuration
        if self.method == "single":
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

            # extract decision variables
            sel = soln.soln_decn[0]

            # shuffle indices for random mating
            self.rng.shuffle(sel)

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
        else:
            raise ValueError("argument 'method' must be either 'single' or 'pareto'")

class SubsetConventionalGenomicSelection(BaseConventionalGenomicSelection):
    """
    Conventional Genomic Selection in a subset search space.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    # use old init

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
        # calculate GEBVs for all individuals
        gebv = gpmod.gebv(gmat).mat

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space = numpy.arange(ntaxa)
        decn_space_lower = numpy.repeat(0, self.nparent)
        decn_space_upper = numpy.repeat(ntaxa-1, self.nparent)

        # construct problem
        prob = GenomicEstimatedBreedingValueSubsetSelectionProblem(
            gebv = gebv,
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

class RealConventionalGenomicSelection(BaseConventionalGenomicSelection):
    """
    Conventional Genomic Selection in a real search space.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    # use old init

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
        # calculate GEBVs for all individuals
        gebv = gpmod.gebv(gmat).mat

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0.0, ntaxa)
        decn_space_upper = numpy.repeat(1.0, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = GenomicEstimatedBreedingValueRealSelectionProblem(
            gebv = gebv,
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
    # use super pareto() function

    ################# Selection Functions ##################
    # use super select() function

class IntegerConventionalGenomicSelection(BaseConventionalGenomicSelection):
    """
    Conventional Genomic Selection in an integer search space.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    # use old init

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
        # calculate GEBVs for all individuals
        gebv = gpmod.gebv(gmat).mat

        # get number of individuals
        ntaxa = gmat.ntaxa

        # get decision space parameters
        decn_space_lower = numpy.repeat(0, ntaxa)
        decn_space_upper = numpy.repeat(ntaxa, ntaxa)
        decn_space = numpy.stack([decn_space_lower,decn_space_upper])

        # construct problem
        prob = GenomicEstimatedBreedingValueIntegerSelectionProblem(
            gebv = gebv,
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
    # use super pareto() function

    ################# Selection Functions ##################
    # use super select() function
