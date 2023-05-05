"""
Module defining selection problems that are by nature binary selection problems.
"""

# list of public objects in this module
__all__ = [
    "BinarySelectionProblemType",
    "check_is_BinarySelectionProblemType"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblemType import SelectionProblemType
from pybrops.opt.prob.BinaryProblemType import BinaryProblemType

# inheritance order not super important here since both abstract
class BinarySelectionProblemType(BinaryProblemType,SelectionProblemType):
    """
    docstring for BinarySelectionProblemType.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for BinarySelectionProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(BinarySelectionProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_BinarySelectionProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type BinarySelectionProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BinarySelectionProblemType):
        raise TypeError("'{0}' must be of type BinarySelectionProblemType.".format(vname))
