"""
Module defining interfaces and error checking routines for subset optimization algorithms.
"""

__all__ = [
    "SubsetOptimizationAlgorithm",
    "check_is_SubsetOptimizationAlgorithm",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Optional
from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
from pybrops.opt.prob.SubsetProblem import SubsetProblem
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class SubsetOptimizationAlgorithm(
        OptimizationAlgorithm,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for optimization algorithms optimizing in subset search spaces.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions in subset search spaces.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def minimize(
            self, 
            prob: SubsetProblem,
            miscout: Optional[dict],
            **kwargs: dict
        ) -> SubsetSolution:
        """
        Minimize an optimization problem.

        Parameters
        ----------
        prob : SubsetProblem
            A subset problem definition object on which to optimize.
        miscout : dict
            Miscellaneous output from the constrained optimizaiont algorithm.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : SubsetSolution
            An object containing the solution to the provided subset problem.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SubsetOptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetOptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,SubsetOptimizationAlgorithm.__name__,type(v).__name__))
