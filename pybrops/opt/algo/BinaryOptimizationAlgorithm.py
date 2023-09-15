"""
Module defining interfaces and error checking routines for binary optimization algorithms.
"""

__all__ = [
    "BinaryOptimizationAlgorithm",
    "check_is_BinaryOptimizationAlgorithm",
]

from abc import ABCMeta, abstractmethod
from typing import Optional
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.opt.prob.BinaryProblem import BinaryProblem
from pybrops.opt.soln.BinarySolution import BinarySolution

class BinaryOptimizationAlgorithm(OptimizationAlgorithm,metaclass=ABCMeta):
    """
    An abstract class for optimization algorithms optimizing in binary search spaces.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions in binary search spaces.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def minimize(
            self, 
            prob: BinaryProblem,
            miscout: Optional[dict],
            **kwargs: dict
        ) -> BinarySolution:
        """
        Minimize an optimization problem.

        Parameters
        ----------
        prob : BinaryProblem
            A binary problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : BinarySolution
            An object containing the solution to the provided binary problem.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_BinaryOptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type BinaryOptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BinaryOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,BinaryOptimizationAlgorithm.__name__,type(v).__name__))
