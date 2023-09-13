Optimization Problems
#####################

Class Family Overview
=====================

All optimization problems in ``PyBrOpS`` are represented using classes derived from the ``Problem`` interface in the ``pybrops.opt.prob`` submodule. All classes derived from the ``Problem`` interface contain metadata pertaining to the decision space, the objective space, and inequality and equality constraints for the optimization. Derivatives of this interface must also implement an evaluation function which is used to evaluate candidate solutions in the decision space. All optimization problems are expressed as minimization problems.

Summary of Optimization Problem Classes
=======================================

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

Deriving Problem Classes
========================

Single objective problem specification
--------------------------------------

.. code-block:: python

    # inherit from RealProblem semi-abstract interface
    class DummySingleObjectiveRealProblem(RealProblem):
        def __init__(
                self, 
                ndecn: Integral, 
                decn_space: Union[numpy.ndarray,None], 
                decn_space_lower: Union[numpy.ndarray,Real,None], 
                decn_space_upper: Union[numpy.ndarray,Real,None], 
                nobj: Integral, 
                obj_wt: Optional[Union[numpy.ndarray,Real]] = None, 
                nineqcv: Optional[Integral] = None, 
                ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None, 
                neqcv: Optional[Integral] = None, 
                eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None, 
                **kwargs: dict
            ):
            """
            Constructor
            """
            super(DummySingleObjectiveRealProblem, self).__init__(
                ndecn, 
                decn_space, 
                decn_space_lower, 
                decn_space_upper, 
                nobj, 
                obj_wt, 
                nineqcv, 
                ineqcv_wt, 
                neqcv, 
                eqcv_wt, 
                **kwargs
            )
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

.. code-block:: python

    class DummyMultiObjectiveRealProblem(RealProblem):
        def __init__(
                self, 
                ndecn: Integral, 
                decn_space: Union[numpy.ndarray,None], 
                decn_space_lower: Union[numpy.ndarray,Real,None], 
                decn_space_upper: Union[numpy.ndarray,Real,None], 
                nobj: Integral, 
                obj_wt: Optional[Union[numpy.ndarray,Real]] = None, 
                nineqcv: Optional[Integral] = None, 
                ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None, 
                neqcv: Optional[Integral] = None, 
                eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None, 
                **kwargs: dict
            ):
            """NA"""
            super(DummyMultiObjectiveRealProblem, self).__init__(
                ndecn, 
                decn_space, 
                decn_space_lower, 
                decn_space_upper, 
                nobj, 
                obj_wt, 
                nineqcv, 
                ineqcv_wt, 
                neqcv, 
                eqcv_wt, 
                **kwargs
            )
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
