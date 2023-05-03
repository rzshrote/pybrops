"""
Module for defining optimization problems with nominal decision variables.
"""

# list of public objects in this module
__all__ = [
    "SubsetProblemType",
    "check_is_SubsetProblemType"
]

# imports
from pybrops.opt.prob.ProblemType import ProblemType

class SubsetProblemType(ProblemType):
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
        super(SubsetProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SubsetProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetProblemType):
        raise TypeError("'{0}' must be of type SubsetProblemType.".format(vname))
