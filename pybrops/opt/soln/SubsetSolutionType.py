"""
Module for defining optimization problem solutions with nominal decision variables.
"""

# list of all public imports in the module
__all__ = [
    "SubsetSolutionType",
    "check_is_SubsetSolutionType"
]

# imports
from pybrops.opt.soln.SolutionType import SolutionType

class SubsetSolutionType(SolutionType):
    """
    Abstract Base Class for all optimization problem solutions with nominal decision variables.
    """
    pass



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SubsetSolutionType(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetSolutionType):
        raise TypeError("'{0}' must be of type SubsetSolution.".format(vname))
