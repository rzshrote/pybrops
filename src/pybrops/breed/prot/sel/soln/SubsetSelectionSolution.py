__all__ = [
    "SubsetSelectionSolution",
]

from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class SubsetSelectionSolution(SubsetSolution,SelectionSolution):
    """
    Class representing subset selection solutions.
    """
    # use implementation from SubsetSolution
    pass



################################## Utilities ###################################
def check_is_SubsetSelectionSolution(v: object, vname: str) -> None:
    """
    Check if object is of type SubsetSelectionSolution, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SubsetSelectionSolution):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,SubsetSelectionSolution.__name__,type(v).__name__))
