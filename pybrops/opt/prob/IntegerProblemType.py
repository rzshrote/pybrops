"""
Module for defining optimization problems with real decision variables.
"""

# list of public objects in this module
__all__ = [
    "IntegerProblemType",
    "check_is_IntegerProblemType"
]

# imports
from pybrops.opt.prob.ProblemType import ProblemType

class IntegerProblemType(ProblemType):
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
        Constructor for IntegerProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(IntegerProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_IntegerProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerProblemType):
        raise TypeError("'{0}' must be of type IntegerProblemType.".format(vname))
