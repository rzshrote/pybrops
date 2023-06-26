__all__ = [
    "IntegerSelectionSolution"
]

from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.opt.soln.IntegerSolution import IntegerSolution

class IntegerSelectionSolution(IntegerSolution,SelectionSolution):
    """
    Class representing subset selection solutions.
    """
    # use implementation from IntegerSolution
    pass



################################## Utilities ###################################
def check_is_IntegerSelectionSolution(v: object, vname: str) -> None:
    """
    Check if object is of type IntegerSelectionSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, IntegerSelectionSolution):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,IntegerSelectionSolution.__name__,type(v).__name__))
