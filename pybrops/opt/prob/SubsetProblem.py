"""
Module for defining optimization problems with nominal decision variables.
"""

# list of public objects in this module
__all__ = [
    "SubsetProblem",
    "check_is_SubsetProblem"
]

# imports
from pybrops.opt.prob.Problem import Problem

class SubsetProblem(Problem):
    """
    Base class for all optimization problems with nominal decision variables.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SubsetProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SubsetProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetProblem):
        raise TypeError("'{0}' must be of type SubsetProblem.".format(vname))
