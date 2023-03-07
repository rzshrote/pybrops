"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "SetSelectionProblem",
    "check_is_SetSelectionProblem"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.opt.prob.SetProblem import SetProblem

# inheritance order not super important here since both abstract
class SetSelectionProblem(SetProblem,SelectionProblem):
    """
    docstring for SetSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SetSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SetSelectionProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_SetSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type SetSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SetSelectionProblem):
        raise TypeError("'{0}' must be of type SetSelectionProblem.".format(vname))
