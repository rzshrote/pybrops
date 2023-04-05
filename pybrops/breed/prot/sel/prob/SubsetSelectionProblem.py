"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "SubsetSelectionProblem",
    "check_is_SubsetSelectionProblem"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.opt.prob.SubsetProblem import SubsetProblem

# inheritance order not super important here since both abstract
class SubsetSelectionProblem(SubsetProblem,SelectionProblem):
    """
    docstring for SubsetSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SubsetSelectionProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SubsetSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetSelectionProblem):
        raise TypeError("'{0}' must be of type SubsetSelectionProblem.".format(vname))
