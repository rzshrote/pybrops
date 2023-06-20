__all__ = [
    "BinarySelectionSolution"
]

from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.opt.soln.BinarySolution import BinarySolution

class BinarySelectionSolution(BinarySolution,SelectionSolution):
    """
    Class representing subset selection solutions.
    """
    # use implementation from BinarySolution
    pass