"""
Module for defining optimization problem solutions with real decision variables.
"""

# list of all public imports in the module
__all__ = [
    "RealSolution",
    "check_is_RealSolution"
]

# imports
from pybrops.opt.soln.Solution import Solution

class RealSolution(Solution):
    """
    Base class for all optimization problem solutions with real decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(RealSolution, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_RealSolution(v: object, vname: str) -> None:
    """
    Check if object is of type RealSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSolution):
        raise TypeError("'{0}' must be of type RealSolution.".format(vname))
