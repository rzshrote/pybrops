"""
Module for defining optimization problems with real decision variables.
"""

# list of public objects in this module
__all__ = [
    "BinaryProblemType",
    "check_is_BinaryProblemType"
]

# imports
from pybrops.opt.prob.ProblemType import ProblemType

class BinaryProblemType(ProblemType):
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
        Constructor for BinaryProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(BinaryProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_BinaryProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type BinaryProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BinaryProblemType):
        raise TypeError("'{0}' must be of type BinaryProblemType.".format(vname))
