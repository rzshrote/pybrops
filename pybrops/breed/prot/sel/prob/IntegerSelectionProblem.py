"""
Module defining selection problems that are by nature set selection problems.
"""

# list of public objects in this module
__all__ = [
    "IntegerSelectionProblem",
    "check_is_IntegerSelectionProblem"
]

# imports
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.opt.prob.IntegerProblemType import IntegerProblemType

# inheritance order not super important here since both abstract
class IntegerSelectionProblem(IntegerProblemType,SelectionProblem):
    """
    docstring for IntegerSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            **kwargs: dict
        ) -> None:
        """
        Constructor for IntegerSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(IntegerSelectionProblem, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_IntegerSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerSelectionProblem):
        raise TypeError("'{0}' must be of type IntegerSelectionProblem.".format(vname))
