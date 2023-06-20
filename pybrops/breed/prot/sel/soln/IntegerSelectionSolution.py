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