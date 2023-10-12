#!/usr/bin/env python3

from numbers import Integral, Real
from typing import Optional, Tuple, Union
import numpy

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

###
### Deriving a Problem class
###

#
# Single objective problem specification
#

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
    # default ``_evaluate`` method inherited from base Problem class

#
# Multi-objective problem specification
#

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
    # default ``_evaluate`` method inherited from base Problem class

###
### Constructing a Problem
###

#
# Construct a single-objective problem
#

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

#
# Construct a multi-objective problem
#

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
### Evaluate candidate solutions
###

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
