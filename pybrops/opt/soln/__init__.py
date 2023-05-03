# order dependent

__all__ = [
    "SolutionType",
    "SubsetSolutionType",
    "RealSolutionType",
    "Solution",
    "SubsetSolution"
]

# base interface
from pybrops.opt.soln import SolutionType

# interfaces derived from Solution
from pybrops.opt.soln import SubsetSolutionType
from pybrops.opt.soln import RealSolutionType

# implementations derived from Solution
from pybrops.opt.soln import Solution

# Implementations derived from SetSolution
from pybrops.opt.soln import SubsetSolution
