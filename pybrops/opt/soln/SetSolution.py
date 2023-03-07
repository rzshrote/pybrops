"""
Module for defining optimization problem solutions with nominal decision variables.
"""

# list of all public imports in the module
__all__ = [
    "SetSolution",
    "check_is_SetSolution"
]

# imports
from pybrops.opt.soln.Solution import Solution

class SetSolution(Solution):
    """
    Base class for all optimization problem solutions with nominal decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SetSolution.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SetSolution, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SetSolution(v: object, vname: str) -> None:
    """
    Check if object is of type SetSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SetSolution):
        raise TypeError("'{0}' must be of type SetSolution.".format(vname))
