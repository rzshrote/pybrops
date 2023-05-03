"""
Module for defining optimization problem solutions with integer decision variables.
"""

# list of all public imports in the module
__all__ = [
    "BinarySolutionType",
    "check_is_BinarySolutionType"
]

# imports
from pybrops.opt.soln.SolutionType import SolutionType

class BinarySolutionType(SolutionType):
    """
    Abstract Base Class for all optimization problem solutions with integer decision variables.
    """
    pass



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_BinarySolutionType(v: object, vname: str) -> None:
    """
    Check if object is of type BinarySolutionType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BinarySolutionType):
        raise TypeError("'{0}' must be of type BinarySolutionType.".format(vname))
