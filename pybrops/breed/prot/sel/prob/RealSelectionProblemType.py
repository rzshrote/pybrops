"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "RealSelectionProblemType",
    "check_is_RealSelectionProblemType"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblemType import SelectionProblemType
from pybrops.opt.prob.RealProblemType import RealProblemType

# inheritance order not super important here since both abstract
class RealSelectionProblemType(RealProblemType,SelectionProblemType):
    """
    docstring for RealSelectionProblemType.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealSelectionProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(RealSelectionProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_RealSelectionProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type RealSelectionProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSelectionProblemType):
        raise TypeError("'{0}' must be of type RealSelectionProblemType.".format(vname))
