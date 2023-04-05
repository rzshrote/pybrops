"""
Module for defining optimization problems with real decision variables.
"""

# list of public objects in this module
__all__ = [
    "RealProblem",
    "check_is_RealProblem"
]

# imports
from pybrops.opt.prob.Problem import Problem

class RealProblem(Problem):
    """
    Base class for all optimization problems with real decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(RealProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_RealProblem(v: object, vname: str) -> None:
    """
    Check if object is of type RealProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealProblem):
        raise TypeError("'{0}' must be of type RealProblem.".format(vname))
