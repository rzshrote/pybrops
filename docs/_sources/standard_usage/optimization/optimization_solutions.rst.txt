Optimization Solutions
######################

Class Family Overview
=====================

All optimization solutions in ``PyBrOpS`` are represented using classes derived from the ``Solution`` interface in the ``pybrops.opt.soln`` submodule. All classes derived from the ``Solution`` interface contain metadata pertaining to the decision space, the objective space, inequality and equality constraints, and the identified solution(s) for the optimization. The purpose of the ``Solution`` family of classes is essentially to store the results of an optimization.

Summary of Optimization Solution Classes
========================================

Optimization solution interfaces and implemented classes can be found in the ``pybrops.opt.soln`` module within PyBrOpS. These classes are summarized below.

.. list-table:: Summary of classes in the ``pybrops.opt.soln`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``Solution``
      - Semi-Abstract
      - Interface for all optimization solution child classes.
    * - ``BinarySolution``
      - Concrete
      - Class representing optimization solutions in binary decision spaces.
    * - ``IntegerSolution``
      - Concrete
      - Class representing optimization solutions in integer decision spaces.
    * - ``RealSolution``
      - Concrete
      - Class representing optimization solutions in real decision spaces.
    * - ``SubsetSolution``
      - Concrete
      - Class representing optimization solutions in subset decision spaces.

Loading Optimization Solution Modules
=====================================

Importing optimization solution classes can be accomplished using the following import statements:

.. code-block:: python

    # import the Solution class (a semi-abstract interface class)
    from pybrops.opt.soln.Solution import Solution

    # import the BinarySolution class (a semi-abstract interface class)
    from pybrops.opt.soln.BinarySolution import BinarySolution

    # import the IntegerSolution class (a semi-abstract interface class)
    from pybrops.opt.soln.IntegerSolution import IntegerSolution

    # import the RealSolution class (a semi-abstract interface class)
    from pybrops.opt.soln.RealSolution import RealSolution

    # import the SubsetSolution class (a semi-abstract interface class)
    from pybrops.opt.soln.SubsetSolution import SubsetSolution

Constructing Solutions
======================

Optimization solutions can be created using their corresponding constructors. The two examples below illustrate the construction of the single- and multi-objective optimization solutions for a real decision space.

Construct a single-objective solution
-------------------------------------

.. code-block:: python

    # solution parameters
    ndecn = 10
    decn_space_lower = numpy.repeat(-1.0, ndecn)
    decn_space_upper = numpy.repeat(1.0, ndecn)
    decn_space = numpy.stack([decn_space_lower, decn_space_upper])
    nobj = 1
    nsoln = 5

    # construct solution
    sosoln = RealSolution(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj,
        obj_wt = 1.0,
        nineqcv = None,
        ineqcv_wt = None,
        neqcv = None,
        eqcv_wt = None,
        nsoln = 5,
        soln_decn = numpy.random.uniform(-1.0, 1.0, (nsoln,ndecn)),
        soln_obj = numpy.random.random((nsoln,nobj)),
        soln_ineqcv = None,
        soln_eqcv = None
    )

Construct a multi-objective solution
------------------------------------

.. code-block:: python

    # solution parameters
    ndecn = 10
    decn_space_lower = numpy.repeat(-1.0, ndecn)
    decn_space_upper = numpy.repeat(1.0, ndecn)
    decn_space = numpy.stack([decn_space_lower, decn_space_upper])
    nobj = 2
    nsoln = 5

    # construct solution
    mosoln = RealSolution(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj,
        obj_wt = 1.0,
        nineqcv = None,
        ineqcv_wt = None,
        neqcv = None,
        eqcv_wt = None,
        nsoln = 5,
        soln_decn = numpy.random.uniform(-1.0, 1.0, (nsoln,ndecn)),
        soln_obj = numpy.random.random((nsoln,nobj)),
        soln_ineqcv = None,
        soln_eqcv = None
    )

Optimization Problem Properties
===============================

All optimization solutions have a set of common properties. The table below summarizes these solution properties.

.. list-table:: Summary of ``Solution`` properties
    :widths: 25 50
    :header-rows: 1

    * - Property
      - Description
    * - ``ndecn``
      - Number of decision variables.
    * - ``decn_space``
      - Decision space boundaries.
    * - ``decn_space_lower``
      - Lower boundary of the decision space.
    * - ``decn_space_upper``
      - Upper boundary of the decision space.
    * - ``nobj``
      - Number of objectives.
    * - ``obj_wt``
      - Objective function weights.
    * - ``nineqcv``
      - Number of inequality constraint violation functions.
    * - ``ineqcv_wt``
      - Inequality constraint violation function weights.
    * - ``neqcv``
      - Number of equality constraint violations.
    * - ``eqcv_wt``
      - Equality constraint violation function weights.
    * - ``nsoln``
      - Number of solutions to the problem.
    * - ``soln_decn``
      - Matrix of solution vectors in the decision space.
    * - ``soln_obj``
      - Solution objective function values.
    * - ``soln_ineqcv``
      - Solution inequality constraint violation function values.
    * - ``soln_eqcv``
      - Solution equality constraint violation function values.

