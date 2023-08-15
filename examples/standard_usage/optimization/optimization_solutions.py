#!/usr/bin/env python3

from numbers import Integral, Real
from typing import Optional, Tuple, Union
import numpy

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

###
### Constructing a Solution
###

#
# Construct a single-objective solution
#

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

#
# Construct a multi-objective solution
#

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
