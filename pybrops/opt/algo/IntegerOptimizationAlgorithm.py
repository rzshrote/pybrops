"""
Module defining interfaces and error checking routines for integer optimization algorithms.
"""

__all__ = [
    "IntegerOptimizationAlgorithm",
    "check_is_IntegerOptimizationAlgorithm",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Optional
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.opt.prob.IntegerProblem import IntegerProblem
from pybrops.opt.soln.IntegerSolution import IntegerSolution

class IntegerOptimizationAlgorithm(
        OptimizationAlgorithm,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for optimization algorithms optimizing in integer search spaces.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions in integer search spaces.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def minimize(
            self, 
            prob: IntegerProblem,
            miscout: Optional[dict],
            **kwargs: dict
        ) -> IntegerSolution:
        """
        Minimize an optimization problem.

        Parameters
        ----------
        prob : IntegerProblem
            A integer problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : IntegerSolution
            An object containing the solution to the provided integer problem.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_IntegerOptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerOptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,IntegerOptimizationAlgorithm.__name__,type(v).__name__))
