# order dependent

__all__ = [
    "Solution",
    "BinarySolution",
    "IntegerSolution",
    "RealSolution",
    "SubsetSolution",
]

# base interface
from pybrops.opt.soln import Solution

# Implementations derived from base interface
from pybrops.opt.soln import BinarySolution
from pybrops.opt.soln import IntegerSolution
from pybrops.opt.soln import RealSolution
from pybrops.opt.soln import SubsetSolution