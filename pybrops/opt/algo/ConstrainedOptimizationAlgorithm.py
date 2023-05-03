"""
Module defining interfaces and associated error checking routines for
constrained optimization algorithms.
"""

from typing import Optional
from pybrops.opt.prob.ProblemType import ProblemType
from pybrops.opt.soln.SolutionType import SolutionType

class ConstrainedOptimizationAlgorithm:
    """
    An abstract class for optimization algorithms.

    The purpose of this abstract class is to provide functionality for:
        1) Optimization of objective functions.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class ConstrainedOptimizationAlgorithm.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(ConstrainedOptimizationAlgorithm, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def minimize(
            self, 
            prob: ProblemType,
            miscout: Optional[dict],
            **kwargs: dict
        ) -> SolutionType:
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



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_ConstrainedOptimizationAlgorithm(v: object, vname: str) -> None:
    """
    Check if object is of type ConstrainedOptimizationAlgorithm. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ConstrainedOptimizationAlgorithm):
        raise TypeError("variable '{0}' must be a ConstrainedOptimizationAlgorithm".format(vname))
