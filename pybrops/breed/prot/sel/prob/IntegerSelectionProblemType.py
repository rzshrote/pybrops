"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "IntegerSelectionProblemType",
    "check_is_IntegerSelectionProblemType"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblemType import SelectionProblemType
from pybrops.opt.prob.IntegerProblemType import IntegerProblemType

# inheritance order not super important here since both abstract
class IntegerSelectionProblemType(IntegerProblemType,SelectionProblemType):
    """
    docstring for IntegerSelectionProblemType.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for IntegerSelectionProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(IntegerSelectionProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_IntegerSelectionProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerSelectionProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerSelectionProblemType):
        raise TypeError("'{0}' must be of type IntegerSelectionProblemType.".format(vname))
