# order dependent

__all__ = [
    "Solution",
    "SubsetSolution",
    "RealSolution",
    "DenseSolution",
    "DenseSubsetSolution"
]

# base interface
from pybrops.opt.soln import Solution

# interfaces derived from Solution
from pybrops.opt.soln import SubsetSolution
from pybrops.opt.soln import RealSolution

# implementations derived from Solution
from pybrops.opt.soln import DenseSolution

# Implementations derived from SetSolution
from pybrops.opt.soln import DenseSubsetSolution
