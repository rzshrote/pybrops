"""
Module containing solution objects
"""

__all__ = [
    "SelectionSolution",
    "MateSelectionSolution",
    "BinarySelectionSolution",
    "IntegerSelectionSolution",
    "RealSelectionSolution",
    "SubsetSelectionSolution",
    "BinaryMateSelectionSolution",
    "IntegerMateSelectionSolution",
    "RealMateSelectionSolution",
    "SubsetMateSelectionSolution"
]

# semi-abstract classes
from pybrops.breed.prot.sel.soln import SelectionSolution
from pybrops.breed.prot.sel.soln import MateSelectionSolution

# concrete classes: tier 1
from pybrops.breed.prot.sel.soln import BinarySelectionSolution
from pybrops.breed.prot.sel.soln import IntegerSelectionSolution
from pybrops.breed.prot.sel.soln import RealSelectionSolution
from pybrops.breed.prot.sel.soln import SubsetSelectionSolution

# concrete classes: tier 2
from pybrops.breed.prot.sel.soln import BinaryMateSelectionSolution
from pybrops.breed.prot.sel.soln import IntegerMateSelectionSolution
from pybrops.breed.prot.sel.soln import RealMateSelectionSolution
from pybrops.breed.prot.sel.soln import SubsetMateSelectionSolution
