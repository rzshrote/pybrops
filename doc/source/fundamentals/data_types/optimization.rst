Optimization Data Types
#######################

One of ``PyBrOpS``'s key features which separates it from other breeding program simulation software is its ability to natively perform single- and multi-objective optimizations. Optimization problems, solutions, and algorithms are encapulated in appropriate object classes. This approach is different from other optimization frameworks which may utilize a function-centric strategy to conduct opitmizations. Infrastructure to define generic opitmization problems, solutions, and algorithms are found in the ``pybrops.opt`` module.

Optimization Problems
*********************

All optimization problems in ``PyBrOpS`` are represented using classes derived from the ``Problem`` interface in the ``pybrops.opt.prob`` submodule. All classes derived from the ``Problem`` interface contain metadata pertaining to the decision space, the objective space, and inequality and equality constraints for the optimization. Derivatives of this interface must also implement an evaluation function which is used to evaluate candidate solutions in the decision space. All optimization problems are expressed as minimization problems.

Optimization Solutions
**********************

All optimization solutions in ``PyBrOpS`` are represented using classes derived from the ``Solution`` interface in the ``pybrops.opt.soln`` submodule. All classes derived from the ``Solution`` interface contain metadata pertaining to the decision space, the objective space, inequality and equality constraints, and the identified solution(s) for the optimization. The purpose of the ``Solution`` family of classes is essentially to store the results of an optimization.

Optimization Algorithms
***********************

All optimization algorithms, regardless of their nature, are represented by classes implementing the ``OptimizationAlgorithm`` interface. Since optimization algorithms are very diverse, the ``OptimizationAlgorithm`` interface only requires a single method to be implemented: the ``minimize`` method.
