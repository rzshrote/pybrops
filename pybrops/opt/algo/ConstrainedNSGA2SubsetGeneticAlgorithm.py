"""
Module implementing an NSGA-II genetic algorithm adapted for subset selection
optimization.
"""

# all public classes and functions available in this module
__all__ = [
    "ConstrainedNSGA2SubsetGeneticAlgorithm"
]

# imports
from numbers import Integral
from numpy.random import Generator,RandomState
from typing import Optional, Union
import numpy
from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.random import global_prng
from pybrops.opt.algo.ConstrainedOptimizationAlgorithm import ConstrainedOptimizationAlgorithm
from pybrops.core.util.pareto import is_pareto_efficient
from pybrops.opt.algo.pymoo_addon import ReducedExchangeCrossover, ReducedExchangeMutation, SubsetRandomSampling
from pybrops.opt.prob.SubsetProblemType import SubsetProblemType
from pybrops.opt.soln.SubsetSolution import SubsetSolution
from pybrops.opt.soln.SubsetSolutionType import SubsetSolutionType
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.termination.max_gen import MaximumGenerationTermination

class ConstrainedNSGA2SubsetGeneticAlgorithm(ConstrainedOptimizationAlgorithm):
    """
    Class implementing an NSGA-II genetic algorithm adapted for subset selection
    optimization. The search space is discrete and nominal in nature.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            ngen: Integral = 250, 
            pop_size: Integral = 100, 
            rng: Union[Generator,RandomState] = global_prng, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for NSGA-II set optimization algorithm.

        Parameters
        ----------
        ngen : int
            Number of generations to evolve population.
        mu : int
            Number of parental candidates to keep in population.
        lamb : int
            Number of progeny to generate.
        M : float
            Length of the chromosome genetic map, in Morgans.
        rng : numpy.random.Generator, numpy.random.RandomState, None
            Random number generator source.
        kwargs : dict
            Additional keyword arguments.
        """
        super(ConstrainedNSGA2SubsetGeneticAlgorithm, self).__init__(**kwargs)
        self.ngen = ngen
        self.pop_size = pop_size
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def ngen(self) -> Integral:
        """Number of generations."""
        return self._ngen
    @ngen.setter
    def ngen(self, value: Integral) -> None:
        """Set number of generations."""
        check_is_Integral(value, "ngen")    # must be int
        check_is_gt(value, "ngen", 0)       # int must be >0
        self._ngen = value
    @ngen.deleter
    def ngen(self) -> None:
        """Delete number of generations."""
        del self._ngen

    @property
    def pop_size(self) -> Integral:
        """Number of individuals in the main chromosome population."""
        return self._pop_size
    @pop_size.setter
    def pop_size(self, value: Integral) -> None:
        """Set number of individuals in the main chromosome population."""
        check_is_Integral(value, "mu")      # must be int
        check_is_gt(value, "mu", 0)         # int must be >0
        self._pop_size = value
    @pop_size.deleter
    def pop_size(self) -> None:
        """Delete number of individuals in the main chromosome population."""
        del self._pop_size

    @property
    def rng(self) -> Union[Generator,RandomState]:
        """Random number generator source."""
        return self._rng
    @rng.setter
    def rng(self, value: Union[Generator,RandomState]) -> None:
        """Set random number generator source."""
        if value is None:
            value = global_prng
        check_is_Generator_or_RandomState(value, "rng")
        self._rng = value
    @rng.deleter
    def rng(self) -> None:
        """Delete random number generator source."""
        del self._rng

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def minimize(
            self, 
            prob: SubsetProblemType,
            miscout: Optional[dict] = None,
            **kwargs: dict
        ) -> SubsetSolutionType:
        """
        Optimize an objective function.

        Parameters
        ----------
        prob : SubsetProblem
            A problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : SubsetSolution
            An object containing the solution to the provided problem.
        """
        # construct the NSGA2 algorithm with custom operators
        algo = NSGA2(
            pop_size = self.pop_size,
            sampling = SubsetRandomSampling(setspace = prob.decn_space),
            crossover = ReducedExchangeCrossover(),
            mutation = ReducedExchangeMutation(setspace = prob.decn_space),
        )

        # optimize the objective function
        res = minimize(
            problem = prob,
            algorithm = algo,
            termination = MaximumGenerationTermination(n_max_gen = self.ngen),
            copy_algorithm = False,
            copy_termination = False
        )

        # extract information from the results object
        if prob.n_obj == 1:
            nsoln = 1
            soln_decn = numpy.stack([res.X])
            soln_obj = numpy.stack([res.F])
            soln_ineqcv = numpy.stack([res.G])
            soln_eqcv = numpy.stack([res.H])
        else:
            nsoln = len(res.X)
            soln_decn = res.X
            soln_obj = res.F
            soln_ineqcv = res.G
            soln_eqcv = res.H

        # construct a solution
        soln = SubsetSolution(
            ndecn = prob.ndecn,
            decn_space = prob.decn_space,
            decn_space_lower = prob.decn_space_lower,
            decn_space_upper = prob.decn_space_upper,
            nobj = prob.nobj,
            obj_wt = prob.obj_wt,
            nineqcv = prob.nineqcv,
            ineqcv_wt = prob.ineqcv_wt,
            neqcv = prob.neqcv,
            eqcv_wt = prob.eqcv_wt,
            nsoln = nsoln,
            soln_decn = soln_decn,
            soln_obj = soln_obj,
            soln_ineqcv = soln_ineqcv,
            soln_eqcv = soln_eqcv
        )

        return soln
