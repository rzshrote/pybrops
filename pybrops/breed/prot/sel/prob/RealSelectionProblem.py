"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "RealSelectionProblem",
    "check_is_RealSelectionProblem"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.opt.prob.RealProblemType import RealProblemType

# inheritance order not super important here since both abstract
class RealSelectionProblem(RealProblemType,SelectionProblem):
    """
    docstring for RealSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(RealSelectionProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_RealSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type RealSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSelectionProblem):
        raise TypeError("'{0}' must be of type RealSelectionProblem.".format(vname))
