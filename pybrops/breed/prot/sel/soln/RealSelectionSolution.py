__all__ = [
    "RealSelectionSolution"
]

from pybrops.breed.prot.sel.soln import SelectionSolution
from pybrops.opt.soln.RealSolution import RealSolution

class RealSelectionSolution(RealSolution,SelectionSolution):
    """
    Class representing subset selection solutions.
    """
    # use implementation from RealSolution
    pass