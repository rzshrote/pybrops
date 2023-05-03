"""
Module for defining optimization problem solutions with real decision variables.
"""

# list of all public imports in the module
__all__ = [
    "RealSolutionType",
    "check_is_RealSolutionType"
]

# imports
from pybrops.opt.soln.SolutionType import SolutionType

class RealSolutionType(SolutionType):
    """
    Abstract Base Class for all optimization problem solutions with real decision variables.
    """
    pass



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_RealSolutionType(v: object, vname: str) -> None:
    """
    Check if object is of type RealSolutionType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSolutionType):
        raise TypeError("'{0}' must be of type RealSolutionType.".format(vname))
