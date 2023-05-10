"""
Module containing optimization algorithms
"""

__all__ = [
    "OptimizationAlgorithm",
    "SteepestAscentSetHillClimber",
    "StochasticAscentSetHillClimber",
    "NSGA2SetGeneticAlgorithm",
    "NSGA3UnityConstraintGeneticAlgorithm"
]

# order dependent imports

# abstract classes
from pybrops.opt.algo import OptimizationAlgorithm

# concrete classes
from pybrops.opt.algo import SteepestAscentSetHillClimber
from pybrops.opt.algo import StochasticAscentSetHillClimber
from pybrops.opt.algo import NSGA2SetGeneticAlgorithm
from pybrops.opt.algo import NSGA3UnityConstraintGeneticAlgorithm
