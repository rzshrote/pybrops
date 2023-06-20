__all__ = [
    "SubsetSelectionSolution"
]

from pybrops.breed.prot.sel.soln.SelectionSolution import SelectionSolution
from pybrops.opt.soln.SubsetSolution import SubsetSolution

class SubsetSelectionSolution(SubsetSolution,SelectionSolution):
    """
    Class representing subset selection solutions.
    """
    # use implementation from SubsetSolution
    pass