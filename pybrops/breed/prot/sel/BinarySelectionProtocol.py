"""
Module defining a general class for binary selection protocols.
"""

from abc import ABCMeta, abstractmethod
from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from numpy.random import Generator, RandomState
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.cfg.BinarySelectionConfiguration import BinarySelectionConfiguration
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.soln.BinarySelectionSolution import BinarySelectionSolution
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.opt.algo.NSGA2BinaryGeneticAlgorithm import NSGA2BinaryGeneticAlgorithm
from pybrops.opt.algo.BinaryGeneticAlgorithm import BinaryGeneticAlgorithm
from pybrops.opt.algo.BinaryOptimizationAlgorithm import BinaryOptimizationAlgorithm, check_is_BinaryOptimizationAlgorithm
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame


class BinarySelectionProtocol(SelectionProtocol,metaclass=ABCMeta):
    """
    Semi-abstract class for creating binary selection protocols.
    """
    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ncross: Integral,
            nparent: Integral,
            nmating: Union[Integral,numpy.ndarray],
            nprogeny: Union[Integral,numpy.ndarray],
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
            soalgo: Optional[OptimizationAlgorithm] = None,
            moalgo: Optional[OptimizationAlgorithm] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedSelectionProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # order dependent assignments
        self.ncross = ncross
        self.nparent = nparent
        self.nmating = nmating
        self.nprogeny = nprogeny
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.obj_trans = obj_trans
        self.obj_trans_kwargs = obj_trans_kwargs
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.ineqcv_trans = ineqcv_trans
        self.ineqcv_trans_kwargs = ineqcv_trans_kwargs
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.eqcv_trans = eqcv_trans
        self.eqcv_trans_kwargs = eqcv_trans_kwargs
        self.ndset_wt = ndset_wt
        self.ndset_trans = ndset_trans
        self.ndset_trans_kwargs = ndset_trans_kwargs
        self.rng = rng
        self.soalgo = soalgo
        self.moalgo = moalgo

    ############################ Object Properties #############################
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
                pop_size = 100  # number of parents in population
            )
            # construct default hillclimber
            # value = SteepestDescentBinaryHillClimber(self.rng)
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
                pop_size = 100  # number of parents in population
            )
        check_is_BinaryOptimizationAlgorithm(value, "moalgo")
        self._moalgo = value

    ############################## Object Methods ##############################

    ########## Optimization Problem Construction ###########
    @abstractmethod
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
        ) -> BinarySelectionProblem:
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
        out : BinarySelectionProblem
            An optimization problem definition.
        """
        raise NotImplementedError("method is abstract")

    ################ Single Objective Solve ################
    def sosolve(
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
        ) -> BinarySelectionSolution:
        """
        Solve the selection problem using a single-objective optimization algorithm.

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
        out : Solution
            A single-objective solution to the posed selection problem.
        """
        # check the number of objectives and raise error if needed
        if self.nobj != 1:
            raise RuntimeError("{0} instance is not single-objective in nature: expected nobj == 1 but received nobj == {1}".format(type(self).__name__,self.nobj))

        # construct the problem
        prob = self.problem(
            pgmat = pgmat, 
            gmat = gmat, 
            ptdf = ptdf, 
            bvmat = bvmat, 
            gpmod = gpmod, 
            t_cur = t_cur, 
            t_max = t_max, 
            **kwargs
        )

        # optimize the problem
        soln = self.soalgo.minimize(
            prob = prob,
            miscout = miscout
        )

        # convert binary solution to binary selection solution
        # this has the exact same metadata as a BinarySolution
        out = BinarySelectionSolution(
            ndecn = soln.ndecn,
            decn_space = soln.decn_space,
            decn_space_lower = soln.decn_space_lower,
            decn_space_upper = soln.decn_space_upper,
            nobj = soln.nobj,
            obj_wt = soln.obj_wt,
            nineqcv = soln.nineqcv,
            ineqcv_wt = soln.ineqcv_wt,
            neqcv = soln.neqcv,
            eqcv_wt = soln.eqcv_wt,
            nsoln = soln.nsoln,
            soln_decn = soln.soln_decn,
            soln_obj = soln.soln_obj,
            soln_ineqcv = soln.soln_ineqcv,
            soln_eqcv = soln.soln_eqcv
        )

        return out

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
        ) -> BinarySelectionSolution:
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
        # check the number of objectives and raise error if needed
        if self.nobj <= 1:
            raise RuntimeError("{0} instance is not multi-objective in nature: expected nobj > 1 but received nobj == {1}".format(type(self).__name__,self.nobj))

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

        # convert binary solution to binary selection solution
        # this has the exact same metadata as a BinarySolution
        out = BinarySelectionSolution(
            ndecn = soln.ndecn,
            decn_space = soln.decn_space,
            decn_space_lower = soln.decn_space_lower,
            decn_space_upper = soln.decn_space_upper,
            nobj = soln.nobj,
            obj_wt = soln.obj_wt,
            nineqcv = soln.nineqcv,
            ineqcv_wt = soln.ineqcv_wt,
            neqcv = soln.neqcv,
            eqcv_wt = soln.eqcv_wt,
            nsoln = soln.nsoln,
            soln_decn = soln.soln_decn,
            soln_obj = soln.soln_obj,
            soln_ineqcv = soln.soln_ineqcv,
            soln_eqcv = soln.soln_eqcv
        )

        return out

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
        ) -> BinarySelectionConfiguration:
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
        out : BinarySelectionConfiguration
            A selection configuration object, requiring all necessary information to mate individuals.
        """
        # if the number of objectives is 1, then we use a single objective algorithm
        if self.nobj == 1:
            # solve the single-objective problem
            sosoln = self.sosolve(
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

            # add solution to miscout if provided.
            if miscout is not None:
                miscout["sosoln"] = sosoln

            # construct a BinarySelectionConfiguration
            selcfg = BinarySelectionConfiguration(
                ncross = self.ncross,
                nparent = self.nparent,
                nmating = self.nmating,
                nprogeny = self.nprogeny,
                pgmat = pgmat,
                xconfig_decn = sosoln.soln_decn[0],
                rng = None
            )

            return selcfg

        # else, we use a multi-objective algorithm and apply a transformation on the non-dominated points to identify a 
        # multi-objective method: objfn_trans returns a multiple values for each
        # selection configuration
        elif self.nobj > 1:
            # solve the single-objective problem
            mosoln = self.mosolve(
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

            # get scores for each of the points along the non-dominated set
            # (nndpts,)
            score = self.ndset_wt * self.ndset_trans(
                mosoln.soln_obj, 
                **self.ndset_trans_kwargs
            )

            # get index of maximum score
            ix = score.argmax()

            # add solution to miscout if provided.
            if miscout is not None:
                miscout["mosoln"] = mosoln

            # construct a BinarySelectionConfiguration
            selcfg = BinarySelectionConfiguration(
                ncross = self.ncross,
                nparent = self.nparent,
                nmating = self.nmating,
                nprogeny = self.nprogeny,
                pgmat = pgmat,
                xconfig_decn = mosoln.soln_decn[ix],
                rng = None
            )

            return selcfg

        # else raise an error as the number of objectives is an illegal value
        else:
            raise ValueError("number of objectives must be greater than zero")
