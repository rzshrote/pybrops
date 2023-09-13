Optimization Algorithms
#######################

Class Family Overview
=====================

All optimization algorithms, regardless of their nature, are represented by classes implementing the ``OptimizationAlgorithm`` interface. Since optimization algorithms are very diverse, the ``OptimizationAlgorithm`` interface only requires a single method to be implemented: the ``minimize`` method.

Summary of Optimization Algorithm Classes
=========================================

.. list-table:: Summary of classes in the ``pybrops.opt.algo`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``OptimizationAlgorithm``
      - Abstract
      - Interface for all optimization algorithm child classes.
    * - ``BinaryOptimizationAlgorithm``
      - Abstract
      - Interface for all optimization algorithm child classes optimizing problems in binary decision spaces.
    * - ``IntegerOptimizationAlgorithm``
      - Abstract
      - Interface for all optimization algorithm child classes optimizing problems in integer decision spaces.
    * - ``RealOptimizationAlgorithm``
      - Abstract
      - Interface for all optimization algorithm child classes optimizing problems in real decision spaces.
    * - ``SubsetOptimizationAlgorithm``
      - Abstract
      - Interface for all optimization algorithm child classes optimizing problems in subset decision spaces.
    * - ``BinaryGeneticAlgorithm``
      - Concrete
      - General genetic algorithm for optimizing single-objective problems in binary decision spaces.
    * - ``IntegerGeneticAlgorithm``
      - Concrete
      - General genetic algorithm for optimizing single-objective problems in integer decision spaces.
    * - ``RealGeneticAlgorithm``
      - Concrete
      - General genetic algorithm for optimizing single-objective problems in real decision spaces.
    * - ``SubsetGeneticAlgorithm``
      - Concrete
      - General genetic algorithm for optimizing single-objective problems in subset decision spaces.
    * - ``NSGA2BinaryGeneticAlgorithm``
      - Concrete
      - NSGA-II genetic algorithm for optimizing multi-objective problems in binary decision spaces.
    * - ``NSGA2IntegerGeneticAlgorithm``
      - Concrete
      - NSGA-II genetic algorithm for optimizing multi-objective problems in integer decision spaces.
    * - ``NSGA2RealGeneticAlgorithm``
      - Concrete
      - NSGA-II genetic algorithm for optimizing multi-objective problems in real decision spaces.
    * - ``NSGA2SubsetGeneticAlgorithm``
      - Concrete
      - NSGA-II genetic algorithm for optimizing multi-objective problems in subset decision spaces.
    * - ``SteepestDescentSubsetHillClimber``
      - Concrete
      - Steepest descent hill-climbing algorithm for optimizing single-objective problems in subset decision spaces.

Loading Optimization Algorithm Modules
======================================

Importing optimization algorithm classes can be accomplished using the following import statements:

.. code-block:: python

    # import OptimizationAlgorithm classes (abstract interface classes)
    from pybrops.opt.algo.OptimizationAlgorithm import OptimizationAlgorithm
    from pybrops.opt.algo.BinaryOptimizationAlgorithm import BinaryOptimizationAlgorithm
    from pybrops.opt.algo.IntegerOptimizationAlgorithm import IntegerOptimizationAlgorithm
    from pybrops.opt.algo.RealOptimizationAlgorithm import RealOptimizationAlgorithm
    from pybrops.opt.algo.SubsetOptimizationAlgorithm import SubsetOptimizationAlgorithm

    # import Genetic Algorithm classes (concrete implementation classes)
    from pybrops.opt.algo.BinaryGeneticAlgorithm import BinaryGeneticAlgorithm
    from pybrops.opt.algo.IntegerGeneticAlgorithm import IntegerGeneticAlgorithm
    from pybrops.opt.algo.RealGeneticAlgorithm import RealGeneticAlgorithm
    from pybrops.opt.algo.SubsetGeneticAlgorithm import SubsetGeneticAlgorithm

    # import NSGA-II Genetic Algorithm classes (concrete implementation classes)
    from pybrops.opt.algo.NSGA2BinaryGeneticAlgorithm import NSGA2BinaryGeneticAlgorithm
    from pybrops.opt.algo.NSGA2IntegerGeneticAlgorithm import NSGA2IntegerGeneticAlgorithm
    from pybrops.opt.algo.NSGA2RealGeneticAlgorithm import NSGA2RealGeneticAlgorithm
    from pybrops.opt.algo.NSGA2SubsetGeneticAlgorithm import NSGA2SubsetGeneticAlgorithm

    # import steepest descent hill climber algorithm (concrete implementation class)
    from pybrops.opt.algo.SteepestDescentSubsetHillClimber import SteepestDescentSubsetHillClimber

Constructing Optimization Algorithms
====================================

Single-objective optimization algorithms
----------------------------------------

.. code-block:: python

    # create a real-valued single-objective genetic algorithm
    soalgo = RealGeneticAlgorithm(ngen = 100, pop_size = 100)

Multi-objective optimization algorithms
---------------------------------------

.. code-block:: python

    # create a real-valued multi-objective genetic algorithm
    moalgo = NSGA2RealGeneticAlgorithm(ngen = 100, pop_size = 100)

Minimizing Optimization Problems
================================

Optimization of single-objective problems
-----------------------------------------

.. code-block:: python

    # minimize single-objective problem and get solution
    sosoln = soalgo.minimize(prob = soprob)


Optimization of multi-objective problems
----------------------------------------

.. code-block:: python

    # minimize multi-objective problem and get solution
    mosoln = moalgo.minimize(prob = moprob)
