"""
Module implementing an NSGA-III genetic algorithm adapted for subset selection
optimization.
"""

# all public classes and functions available in this module
__all__ = [
    "NSGA3SubsetGeneticAlgorithm",
]

# imports
from numbers import Integral
from typing import Optional
from typing import Union
import numpy
from numpy.random import Generator
from numpy.random import RandomState
from pymoo.algorithms.moo.nsga3 import NSGA3
from pymoo.optimize import minimize
from pymoo.termination.max_gen import MaximumGenerationTermination
from pymoo.util.ref_dirs import get_reference_directions

from pybrops.core.error.error_type_numpy import check_is_Generator_or_RandomState
from pybrops.core.error.error_type_python import check_is_Integral
from pybrops.core.error.error_type_python import check_is_dict
from pybrops.core.error.error_value_python import check_is_gt
from pybrops.core.random.prng import global_prng
from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm
from pybrops.opt.algo.pymoo_addon import ReducedExchangeCrossover
from pybrops.opt.algo.pymoo_addon import ReducedExchangeMutation
from pybrops.opt.algo.pymoo_addon import SubsetRandomSampling
from pybrops.opt.prob.SubsetProblem import SubsetProblem
from pybrops.opt.prob.SubsetProblem import check_SubsetProblem_is_multi_objective
from pybrops.opt.prob.SubsetProblem import check_is_SubsetProblem
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class NSGA3SubsetGeneticAlgorithm(SubsetOptimizationAlgorithm):
    """
    Class implementing an NSGA-III genetic algorithm adapted for subset selection
    optimization. The search space is discrete and nominal in nature.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            ngen: Integral = 250, 
            pop_size: Integral = 100, 
            nrefpts: Optional[Integral] = None, 
            rng: Union[Generator,RandomState] = global_prng, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for NSGA-III subset optimization algorithm.

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
        self.ngen = ngen
        self.pop_size = pop_size
        self.nrefpts = nrefpts
        self.rng = rng

    ############################ Object Properties #############################
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

    @property
    def nrefpts(self) -> Integral:
        """Number of reference points to be used for optimization."""
        return self._nrefpts
    @nrefpts.setter
    def nrefpts(self, value: Integral) -> None:
        """Set the number of reference points to be used for optimization."""
        # if None, default to pop_size - 4 to allow freedom of movement.
        if value is None:
            value = self.pop_size
        check_is_Integral(value, "nrefpts") # must be int
        check_is_gt(value, "nrefpts", 0)    # int must be >0
        self._nrefpts = value

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

    ############################## Object Methods ##############################
    def minimize(
            self, 
            prob: SubsetProblem,
            miscout: Optional[dict] = None,
            **kwargs: dict
        ) -> SubsetSolution:
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
        # type checks
        check_is_SubsetProblem(prob, "prob")
        check_SubsetProblem_is_multi_objective(prob, "prob")
        if miscout is not None:
            check_is_dict(miscout, "miscout")
        
        # create the reference directions to be used for the optimization
        ref_dirs = get_reference_directions(
            "das-dennis", 
            n_dim = prob.nobj, 
            n_points = self.nrefpts,
        )

        # construct the NSGA2 algorithm with custom operators
        algo = NSGA3(
            ref_dirs = ref_dirs,
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
