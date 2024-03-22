__all__ = [
    "RealSelectionSolution",
]

from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.opt.soln.RealSolution import RealSolution

class RealSelectionSolution(
        RealSolution,
        SelectionSolution,
    ):
    """
    Class representing subset selection solutions.
    """
    # use implementation from RealSolution
    pass



################################## Utilities ###################################
def check_is_RealSelectionSolution(v: object, vname: str) -> None:
    """
    Check if object is of type RealSelectionSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, RealSelectionSolution):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,RealSelectionSolution.__name__,type(v).__name__))
