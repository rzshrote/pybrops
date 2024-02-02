"""
Module defining interfaces and associated error checking routines for
constrained optimization algorithms.
"""

__all__ = [
    "OptimizationAlgorithm",
    "check_is_OptimizationAlgorithm",
]

from abc import ABCMeta, abstractmethod
from typing import Optional
from pybrops.opt.prob.Problem import Problem
from pybrops.opt.soln.Solution import Solution

class OptimizationAlgorithm(metaclass=ABCMeta):
    """
    An abstract class for optimization algorithms.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def minimize(
            self, 
            prob: Problem,
            miscout: Optional[dict],
            **kwargs: dict
        ) -> Solution:
        """
        Minimize an optimization problem.

        Parameters
        ----------
        prob : Problem
            A problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Solution
            An object containing the solution to the provided problem.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_OptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type OptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, OptimizationAlgorithm):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,OptimizationAlgorithm.__name__,type(v).__name__))
