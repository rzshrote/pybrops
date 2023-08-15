"""
Module containing optimization algorithms
"""

__all__ = [
    "pymoo_addon",
    "OptimizationAlgorithm",
    "BinaryOptimizationAlgorithm",
    "IntegerOptimizationAlgorithm",
    "RealOptimizationAlgorithm",
    "SubsetOptimizationAlgorithm",
    "BinaryGeneticAlgorithm",
    "IntegerGeneticAlgorithm",
    "RealGeneticAlgorithm",
    "SubsetGeneticAlgorithm",
    "NSGA2BinaryGeneticAlgorithm",
    "NSGA2IntegerGeneticAlgorithm",
    "NSGA2RealGeneticAlgorithm",
    "NSGA2SubsetGeneticAlgorithm",
    "SteepestDescentSubsetHillClimber",

    "UnconstrainedOptimizationAlgorithm",
    "UnconstrainedSteepestAscentSetHillClimber",
    "UnconstrainedStochasticAscentSetHillClimber",
    "UnconstrainedNSGA2SetGeneticAlgorithm",
    "UnconstrainedNSGA3UnityConstraintGeneticAlgorithm"
]

# order dependent imports

# utility imports (tier 0)
from pybrops.opt.algo import pymoo_addon

# abstract classes (tier 1)
from pybrops.opt.algo import OptimizationAlgorithm

# abstract classes (tier 2)
from pybrops.opt.algo import BinaryOptimizationAlgorithm
from pybrops.opt.algo import IntegerOptimizationAlgorithm
from pybrops.opt.algo import RealOptimizationAlgorithm
from pybrops.opt.algo import SubsetOptimizationAlgorithm

# concrete classes (tier 3)
# GAs
from pybrops.opt.algo import BinaryGeneticAlgorithm
from pybrops.opt.algo import IntegerGeneticAlgorithm
from pybrops.opt.algo import RealGeneticAlgorithm
from pybrops.opt.algo import SubsetGeneticAlgorithm
from pybrops.opt.algo import NSGA2BinaryGeneticAlgorithm
from pybrops.opt.algo import NSGA2IntegerGeneticAlgorithm
from pybrops.opt.algo import NSGA2RealGeneticAlgorithm
from pybrops.opt.algo import NSGA2SubsetGeneticAlgorithm
# other
from pybrops.opt.algo import SteepestDescentSubsetHillClimber


# old classes
from pybrops.opt.algo import UnconstrainedOptimizationAlgorithm
from pybrops.opt.algo import UnconstrainedSteepestAscentSetHillClimber
from pybrops.opt.algo import UnconstrainedStochasticAscentSetHillClimber
from pybrops.opt.algo import UnconstrainedNSGA2SetGeneticAlgorithm
from pybrops.opt.algo import UnconstrainedNSGA3UnityConstraintGeneticAlgorithm
