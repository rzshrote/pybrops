"""
Module defining interfaces and error checking routines for real optimization algorithms.
"""

__all__ = [
    "RealOptimizationAlgorithm",
    "check_is_RealOptimizationAlgorithm",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Optional
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.opt.prob.RealProblem import RealProblem
from pybrops.opt.soln.RealSolution import RealSolution

class RealOptimizationAlgorithm(
        OptimizationAlgorithm,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for optimization algorithms optimizing in real search spaces.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions in real search spaces.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def minimize(
            self, 
            prob: RealProblem,
            miscout: Optional[dict],
            **kwargs: dict
        ) -> RealSolution:
        """
        Minimize an optimization problem.

        Parameters
        ----------
        prob : RealProblem
            A real problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : RealSolution
            An object containing the solution to the provided real problem.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_RealOptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type RealOptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,RealOptimizationAlgorithm.__name__,type(v).__name__))
