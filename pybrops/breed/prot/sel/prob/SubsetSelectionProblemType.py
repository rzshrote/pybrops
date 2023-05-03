"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "SubsetSelectionProblemType",
    "check_is_SubsetSelectionProblemType"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblemType import SelectionProblemType
from pybrops.opt.prob.SubsetProblemType import SubsetProblemType

# inheritance order not super important here since both abstract
class SubsetSelectionProblemType(SubsetProblemType,SelectionProblemType):
    """
    docstring for SubsetSelectionProblemType.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetSelectionProblemType.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SubsetSelectionProblemType, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SubsetSelectionProblemType(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetSelectionProblemType, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetSelectionProblemType):
        raise TypeError("'{0}' must be of type SubsetSelectionProblemType.".format(vname))
