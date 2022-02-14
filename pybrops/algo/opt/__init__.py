"""
Module containing optimization algorithms
"""
# order dependent imports

# abstract classes
from . import OptimizationAlgorithm

# concrete classes
from . import SteepestAscentSetHillClimber
from . import StochasticAscentSetHillClimber
from . import NSGA2SetGeneticAlgorithm
