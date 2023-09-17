#!/usr/bin/env python3

from numbers import Integral, Real
from typing import Optional, Tuple, Union
import numpy
from pybrops.opt.prob.RealProblem import RealProblem

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

#
# Create dummy single- and multi-objective problem specifications
#

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

#
# Constructing test problems
#

# Construct a single-objective problem
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

# Construct a multi-objective problem
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

###
### Constructing an optimization algorithm
###

#
# Single-objective optimization algorithm construction
#

# create a real-valued single-objective genetic algorithm
soalgo = RealGeneticAlgorithm(ngen = 100, pop_size = 100)

#
# Multi-objective optimization algorithm construction
#

# create a real-valued multi-objective genetic algorithm
moalgo = NSGA2RealGeneticAlgorithm(ngen = 100, pop_size = 100)

###
### Optimizating Problems
###

#
# Single-objective optimization
#

# minimize single-objective problem and get solution
sosoln = soalgo.minimize(prob = soprob)

#
# Multi-objective optimzation
#

# minimize multi-objective problem and get solution
mosoln = moalgo.minimize(prob = moprob)
