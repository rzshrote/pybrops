Optimization Solutions
######################

Class Family Overview
=====================

All optimization solutions in ``PyBrOpS`` are represented using classes derived from the ``Solution`` interface in the ``pybrops.opt.soln`` submodule. All classes derived from the ``Solution`` interface contain metadata pertaining to the decision space, the objective space, inequality and equality constraints, and the identified solution(s) for the optimization. The purpose of the ``Solution`` family of classes is essentially to store the results of an optimization.

Summary of Optimization Solution Classes
========================================

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
