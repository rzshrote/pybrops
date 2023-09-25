Optimization Problems
#####################

Class Family Overview
=====================

All optimization problems in ``PyBrOpS`` are represented using classes derived from the ``Problem`` interface in the ``pybrops.opt.prob`` submodule. All classes derived from the ``Problem`` interface contain metadata pertaining to the decision space, the objective space, and inequality and equality constraints for the optimization. Derivatives of this interface must also implement an evaluation function which is used to evaluate candidate solutions in the decision space. All optimization problems are expressed as minimization problems.

Summary of Optimization Problem Classes
=======================================

Definitions of optimization problem interfaces within PyBrOpS can be found in the ``pybrops.opt.prob`` module. Within this module are several semi-abstract class definitions which are summarized below.

.. list-table:: Summary of classes in the ``pybrops.opt.prob`` module
    :widths: 25 20 50
    :header-rows: 1

    * - Class Name
      - Class Type
      - Class Description
    * - ``Problem``
      - Semi-Abstract
      - Interface for all optimization problem child classes.
    * - ``BinaryProblem``
      - Semi-Abstract
      - Interface for all optimization problem child classes representing problems in binary decision spaces.
    * - ``IntegerProblem``
      - Semi-Abstract
      - Interface for all optimization problem child classes representing problems in integer decision spaces.
    * - ``RealProblem``
      - Semi-Abstract
      - Interface for all optimization problem child classes representing problems in real decision spaces.
    * - ``SubsetProblem``
      - Semi-Abstract
      - Interface for all optimization problem child classes representing problems in subset decision spaces.

Loading Optimization Problem Modules
====================================

Importing optimization problem classes can be accomplished using the following import statements:

.. code-block:: python

    # import the Problem class (a semi-abstract interface class)
    from pybrops.opt.prob.Problem import Problem

    # import the BinaryProblem class (a semi-abstract interface class)
    from pybrops.opt.prob.BinaryProblem import BinaryProblem

    # import the IntegerProblem class (a semi-abstract interface class)
    from pybrops.opt.prob.IntegerProblem import IntegerProblem

    # import the RealProblem class (a semi-abstract interface class)
    from pybrops.opt.prob.RealProblem import RealProblem

    # import the SubsetProblem class (a semi-abstract interface class)
    from pybrops.opt.prob.SubsetProblem import SubsetProblem

Optimization Problem Properties
===============================

PyMOO specific properties
-------------------------

All optimization problem definitions inherit from PyMOO's Problem definition. For information regarding these properties, see `PyMOO's Problem Definition <https://pymoo.org/problems/definition.html>`_ webpage.

PyBrOpS specific properties
---------------------------

To create a more stable interface, PyBrOpS defines several properties within its Problem definition. These additional properties are summarized below.

.. list-table:: Summary of PyBrOpS specific ``Problem`` properties
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

Deriving Problem Classes
========================

If one desires to create a new optimization problem, one must define a custom Problem class which implements two abstract methods: ``evalfn`` and ``_evaluate``. The ``evalfn`` method is defined and utilized by PyBrOpS internals, while the ``_evaluate`` method is utilized by PyMOO internals. The PyMOO method has a more permissive set of parameters whereas the PyBrOpS method adheres to a more strict interface.

The following two examples illustrate the definition of single- and multi-objective optimization problems.

Single objective problem specification
--------------------------------------

The class definition below defines an optimization problem for the sphere problem. The sphere problem is defined as:

.. math::

    \arg \min_{\mathbf{x}} f(\mathbf{x}) = \sum_{i=1}^{n_{decn}} x_{i}^{2}

.. code-block:: python

    # inherit from RealProblem semi-abstract interface
    class DummySingleObjectiveRealProblem(RealProblem):
        # inherit init since it is unchanged
        ### method required by PyBrOpS interface ###
        def evalfn(
                self, 
                x: numpy.ndarray, 
                *args: tuple, 
                **kwargs: dict
            ) -> Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
            """
            Evaluate a candidate solution for the "Sphere Problem".

            This calculates three vectors which are to be minimized:

            .. math::

                \\mathbf{v_{obj}} = \\mathbf{w_{obj} \\odot F_{obj}(x)} \\
                \\mathbf{v_{ineqcv}} = \\mathbf{w_{ineqcv} \\odot G_{ineqcv}(x)} \\
                \\mathbf{v_{eqcv}} = \\mathbf{w_{eqcv} \\odot H_{eqcv}(x)}
            
            Parameters
            ----------
            x : numpy.ndarray
                A candidate solution vector of shape ``(ndecn,)``.
            args : tuple
                Additional non-keyword arguments.
            kwargs : dict
                Additional keyword arguments.
            
            Returns
            -------
            out : tuple
                A tuple ``(obj, ineqcv, eqcv)``.
                
                Where:
                
                - ``obj`` is a numpy.ndarray of shape ``(nobj,)`` that contains 
                    objective function evaluations.
                - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                    inequality constraint violation values.
                - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                    equality constraint violation values.
            """
            obj = self.obj_wt * numpy.sum(x**2, axis=0, keepdims=True)
            ineqcv = self.ineqcv_wt * numpy.zeros(self.nineqcv)
            eqcv = self.eqcv_wt * numpy.zeros(self.neqcv)
            return obj, ineqcv, eqcv
        ### method required by PyMOO interface ###
        def _evaluate(
                self, 
                x: numpy.ndarray, 
                out: dict, 
                *args: tuple, 
                **kwargs: dict
            ) -> None:
            """
            Evaluate a set of candidate solutions for the "Sphere Problem".

            Parameters
            ----------
            x : numpy.ndarray
                A candidate solution vector of shape ``(nsoln,ndecn)``.
                Where ``nsoln`` is the number of candidates solutions and ``ndecn``
                is the number of decision variables.
            out : dict
                Dictionary to which to output function evaluations.
            args : tuple
                Additional arguments.
            kwargs : dict
                Additional keyword arguments.
            """
            # the PyMOO interface demands acceptance of 1d or 2d numpy.ndarrays
            # this handles either case
            if x.ndim == 1:
                # get evaluations
                vals = self.evalfn(x, *args, **kwargs)
                # create temporary dictionary
                tmp = {key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0}
                # update output dictionary
                out.update(tmp)
            else:
                # create lists for accumulating variables
                objs = []
                ineqcvs = []
                eqcvs = []
                # for each row in x
                for v in x:
                    # get evaluations
                    obj, ineqcv, eqcv = self.evalfn(v, *args, **kwargs)
                    # append values to lists
                    objs.append(obj)
                    ineqcvs.append(ineqcv)
                    eqcvs.append(eqcv)
                # stack outputs
                objs = numpy.stack(objs)
                ineqcvs = numpy.stack(ineqcvs)
                eqcvs = numpy.stack(eqcvs)
                # create temporary dictionary
                tmp = {key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0}
                # update output dictionary
                out.update(tmp)


Multi-objective problem specification
-------------------------------------

The class definition below defines an optimization problem for a multi-objective, dual sphere problem. The dual sphere problem is defined as:

.. math::

    \arg \min_{\mathbf{x}} f_{1}(\mathbf{x}) = \sum_{i=1}^{n_{decn}} x_{i}^{2}
    
    \arg \min_{\mathbf{x}} f_{2}(\mathbf{x}) = \sum_{i=1}^{n_{decn}} (1-x_{i})^{2}

.. code-block:: python

    class DummyMultiObjectiveRealProblem(RealProblem):
        # inherit init since it is unchanged
        ### method required by PyBrOpS interface ###
        def evalfn(
                self, 
                x: numpy.ndarray, 
                *args: tuple, 
                **kwargs: dict
            ) -> Tuple[numpy.ndarray,numpy.ndarray,numpy.ndarray]:
            """
            Evaluate a candidate solution for a Dual Sphere Problem.

            This calculates three vectors which are to be minimized:

            .. math::

                \\mathbf{v_{obj}} = \\mathbf{w_{obj} \\odot F_{obj}(x)} \\
                \\mathbf{v_{ineqcv}} = \\mathbf{w_{ineqcv} \\odot G_{ineqcv}(x)} \\
                \\mathbf{v_{eqcv}} = \\mathbf{w_{eqcv} \\odot H_{eqcv}(x)}
            
            Parameters
            ----------
            x : numpy.ndarray
                A candidate solution vector of shape ``(ndecn,)``.
            args : tuple
                Additional non-keyword arguments.
            kwargs : dict
                Additional keyword arguments.
            
            Returns
            -------
            out : tuple
                A tuple ``(obj, ineqcv, eqcv)``.
                
                Where:
                
                - ``obj`` is a numpy.ndarray of shape ``(nobj,)`` that contains 
                    objective function evaluations.
                - ``ineqcv`` is a numpy.ndarray of shape ``(nineqcv,)`` that contains 
                    inequality constraint violation values.
                - ``eqcv`` is a numpy.ndarray of shape ``(neqcv,)`` that contains 
                    equality constraint violation values.
            """
            obj = self.obj_wt * numpy.array([numpy.sum(x**2), numpy.sum((1-x)**2)], dtype=float)
            ineqcv = self.ineqcv_wt * numpy.zeros(self.nineqcv)
            eqcv = self.eqcv_wt * numpy.zeros(self.neqcv)
            return obj, ineqcv, eqcv
        ### method required by PyMOO interface ###
        def _evaluate(
                self, 
                x: numpy.ndarray, 
                out: dict, 
                *args: tuple, 
                **kwargs: dict
            ) -> None:
            """
            Evaluate a set of candidate solutions for a Dual Sphere Problem.

            Parameters
            ----------
            x : numpy.ndarray
                A candidate solution vector of shape ``(nsoln,ndecn)``.
                Where ``nsoln`` is the number of candidates solutions and ``ndecn``
                is the number of decision variables.
            out : dict
                Dictionary to which to output function evaluations.
            args : tuple
                Additional arguments.
            kwargs : dict
                Additional keyword arguments.
            """
            if x.ndim == 1:
                # get evaluations
                vals = self.evalfn(x, *args, **kwargs)
                # create temporary dictionary
                tmp = {key:val for key,val in zip(["F","G","H"],vals) if len(val) > 0}
                # update output dictionary
                out.update(tmp)
            else:
                # create lists for accumulating variables
                objs = []
                ineqcvs = []
                eqcvs = []
                # for each row in x
                for v in x:
                    # get evaluations
                    obj, ineqcv, eqcv = self.evalfn(v, *args, **kwargs)
                    # append values to lists
                    objs.append(obj)
                    ineqcvs.append(ineqcv)
                    eqcvs.append(eqcv)
                # stack outputs
                objs = numpy.stack(objs)
                ineqcvs = numpy.stack(ineqcvs)
                eqcvs = numpy.stack(eqcvs)
                # create temporary dictionary
                tmp = {key:val for key,val in zip(["F","G","H"],[obj,ineqcv,eqcv]) if val.shape[1] > 0}
                # update output dictionary
                out.update(tmp)

Constructing Problems
=====================

Optimization problems can be created using their corresponding constructors. The two examples below illustrate the construction of the single- and multi-objective optimization problems defined above.

Construct a single-objective problem
------------------------------------

.. code-block:: python

    # problem parameters
    ndecn = 10
    decn_space_lower = numpy.repeat(-1.0, ndecn)
    decn_space_upper = numpy.repeat(1.0, ndecn)
    decn_space = numpy.stack([decn_space_lower, decn_space_upper])
    nobj = 1

    soprob = DummySingleObjectiveRealProblem(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj
    )

Construct a multi-objective problem
-----------------------------------

.. code-block:: python

    # problem parameters
    ndecn = 10
    decn_space_lower = numpy.repeat(-1.0, ndecn)
    decn_space_upper = numpy.repeat(1.0, ndecn)
    decn_space = numpy.stack([decn_space_lower, decn_space_upper])
    nobj = 2

    moprob = DummyMultiObjectiveRealProblem(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj
    )

Evaluating candidate solutions
==============================

Single, candidate solutions can be evaluated using the ``evalfn`` method. In ``PyBrOpS``, candidate solution vectors must be evaluated individually, as is defined by the ``Problem`` interface.  

.. code-block:: python

    # create a random vector in the decision space
    cand = numpy.random.uniform(-1.0, 1.0, ndecn)

    # evaluate the candidate solution
    obj1, ineqcv1, eqcv1 = soprob.evalfn(cand)  # obj1 is length 1
    obj2, ineqcv2, eqcv2 = moprob.evalfn(cand)  # obj2 is length 2

    # create 5 random vectors in the decision space
    cands = numpy.random.uniform(-1.0, 1.0, (5,ndecn))

    # evaluate the candidate solutions and unpack into components
    evaluations = [soprob.evalfn(cand) for cand in cands]
    obj1, ineqcv1, eqcv1 = zip(*evaluations)
    obj1 = numpy.stack(obj1)        # convert list to 2d array
    ineqcv1 = numpy.stack(ineqcv1)  # convert list to 2d array
    eqcv1 = numpy.stack(eqcv1)      # convert list to 2d array

    evaluations = [moprob.evalfn(cand) for cand in cands]
    obj2, ineqcv2, eqcv2 = zip(*evaluations)
    obj2 = numpy.stack(obj2)        # convert list to 2d array
    ineqcv2 = numpy.stack(ineqcv2)  # convert list to 2d array
    eqcv2 = numpy.stack(eqcv2)      # convert list to 2d array
