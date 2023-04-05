from numbers import Integral, Number
from typing import Callable, Iterable, Sequence
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.DenseSolution import DenseSolution, check_is_DenseSolution

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def ndecn():
    yield 4

@pytest.fixture
def decn_space_lower():
    yield numpy.array([1,2,3,4], dtype=float)

@pytest.fixture
def decn_space_upper():
    yield numpy.array([5,6,7,8], dtype=float)

@pytest.fixture
def decn_space(decn_space_lower, decn_space_upper):
    yield numpy.stack([decn_space_lower, decn_space_upper])

@pytest.fixture
def nobj():
    yield 2

@pytest.fixture
def obj_wt():
    yield numpy.array([1,-1], dtype=float)

@pytest.fixture
def nineqcv():
    yield 0

@pytest.fixture
def ineqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def neqcv():
    yield 0

@pytest.fixture
def eqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def nsoln():
    yield 5

@pytest.fixture
def soln_decn(nsoln, ndecn, decn_space_lower, decn_space_upper):
    yield numpy.random.randint(decn_space_lower, decn_space_upper, (nsoln,ndecn))

@pytest.fixture
def soln_obj(nsoln, nobj):
    yield numpy.random.random((nsoln,nobj))

@pytest.fixture
def soln_ineqcv(nsoln, nineqcv):
    yield numpy.random.random((nsoln,nineqcv))

@pytest.fixture
def soln_eqcv(nsoln, neqcv):
    yield numpy.random.random((nsoln,neqcv))

@pytest.fixture
def prob(
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
    nsoln,
    soln_decn,
    soln_obj,
    soln_ineqcv,
    soln_eqcv
):
    yield DenseSolution(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj,
        obj_wt = obj_wt,
        nineqcv = nineqcv,
        ineqcv_wt = ineqcv_wt,
        neqcv = neqcv,
        eqcv_wt = eqcv_wt,
        nsoln = nsoln,
        soln_decn = soln_decn,
        soln_obj = soln_obj,
        soln_ineqcv = soln_ineqcv,
        soln_eqcv = soln_eqcv
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_DenseSolution_docstring():
    assert_docstring(DenseSolution)

################################################################################
########################### Test concrete properties ###########################
################################################################################

#############
### ndecn ###
#############
def test_ndecn_is_concrete():
    assert_concrete_property(DenseSolution, "ndecn")

def test_ndecn_fget(prob):
    assert isinstance(prob.ndecn, Integral)

def test_ndecn_fset(prob):
    with not_raises(TypeError):
        prob.ndecn = numpy.int8(1)
    with not_raises(TypeError):
        prob.ndecn = numpy.int16(1)
    with not_raises(TypeError):
        prob.ndecn = numpy.int32(1)
    with not_raises(TypeError):
        prob.ndecn = numpy.int64(1)
    with not_raises(TypeError):
        prob.ndecn = int(1)

def test_ndecn_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ndecn = None
    with pytest.raises(TypeError):
        prob.ndecn = "string"
    with pytest.raises(TypeError):
        prob.ndecn = float(1.0)

def test_ndecn_fset_ValueError(prob):
    with pytest.raises(ValueError):
        prob.ndecn = int(0)
    with pytest.raises(ValueError):
        prob.ndecn = int(-1)

def test_ndecn_fdel(prob):
    del prob.ndecn
    with pytest.raises(AttributeError):
        prob.ndecn

##################
### decn_space ###
##################
def test_decn_space_is_concrete():
    assert_concrete_property(DenseSolution, "decn_space")

def test_decn_space_fget(prob, ndecn):
    assert isinstance(prob.decn_space, numpy.ndarray) or (prob.decn_space is None)
    if prob.decn_space is not None:
        assert prob.decn_space.shape == (2,ndecn)

def test_decn_space_fset(prob, decn_space_lower, decn_space_upper):
    with not_raises(TypeError):
        prob.decn_space = numpy.stack([decn_space_lower,decn_space_upper])
    with not_raises(TypeError):
        prob.decn_space = None

def test_decn_space_fset_TypeError(prob, ndecn):
    with pytest.raises(TypeError):
        prob.decn_space = "string"
    with pytest.raises(TypeError):
        prob.decn_space = int(1)
    with pytest.raises(TypeError):
        prob.decn_space = ndecn * [1]

def test_decn_space_fset_ValueError(prob, decn_space_lower, decn_space_upper):
    with pytest.raises(ValueError):
        prob.decn_space = numpy.array([1])
    with pytest.raises(ValueError):
        prob.decn_space = decn_space_lower
    with pytest.raises(ValueError):
        prob.decn_space = numpy.concatenate([decn_space_lower,decn_space_upper])
    with pytest.raises(ValueError):
        prob.decn_space = numpy.stack([
            decn_space_lower,decn_space_lower,
            decn_space_upper,decn_space_upper
        ])
    with pytest.raises(ValueError):
        prob.decn_space = numpy.stack([
            numpy.concatenate([decn_space_lower,decn_space_lower]),
            numpy.concatenate([decn_space_upper,decn_space_upper])
        ])

def test_decn_space_fdel(prob):
    del prob.decn_space
    with pytest.raises(AttributeError):
        prob.decn_space

########################
### decn_space_lower ###
########################
def test_decn_space_lower_is_concrete():
    assert_concrete_property(DenseSolution, "decn_space_lower")

def test_decn_space_lower_fget(prob):
    assert isinstance(prob.decn_space_lower, numpy.ndarray) or (prob.decn_space_lower is None)
    assert len(prob.decn_space_lower) == prob.ndecn

def test_decn_space_lower_fset(prob, decn_space_lower):
    with not_raises(TypeError):
        prob.decn_space_lower = decn_space_lower
    with not_raises(TypeError):
        prob.decn_space_lower = decn_space_lower[0]
    with not_raises(TypeError):
        prob.decn_space_lower = None

def test_decn_space_lower_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.decn_space_lower = "string"

def test_decn_space_lower_fset_ValueError(prob, decn_space_lower):
    with pytest.raises(ValueError):
        prob.decn_space_lower = numpy.concatenate([decn_space_lower, decn_space_lower])

def test_decn_space_lower_fdel(prob):
    del prob.decn_space_lower
    with pytest.raises(AttributeError):
        prob.decn_space_lower

########################
### decn_space_upper ###
########################
def test_decn_space_upper_is_concrete():
    assert_concrete_property(DenseSolution, "decn_space_upper")

def test_decn_space_upper_fget(prob):
    assert isinstance(prob.decn_space_upper, numpy.ndarray) or (prob.decn_space_upper is None)
    assert len(prob.decn_space_upper) == prob.ndecn

def test_decn_space_upper_fset(prob, decn_space_upper):
    with not_raises(TypeError):
        prob.decn_space_upper = decn_space_upper
    with not_raises(TypeError):
        prob.decn_space_upper = decn_space_upper[0]
    with not_raises(TypeError):
        prob.decn_space_upper = None

def test_decn_space_upper_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.decn_space_upper = "string"

def test_decn_space_upper_fset_ValueError(prob, decn_space_upper):
    with pytest.raises(ValueError):
        prob.decn_space_upper = numpy.concatenate([decn_space_upper, decn_space_upper])

def test_decn_space_upper_fdel(prob):
    del prob.decn_space_upper
    with pytest.raises(AttributeError):
        prob.decn_space_upper

############
### nobj ###
############
def test_nobj_is_concrete():
    assert_concrete_property(DenseSolution, "nobj")

def test_nobj_fget(prob):
    assert isinstance(prob.nobj, Integral)

def test_nobj_fset(prob):
    with not_raises(TypeError):
        prob.nobj = numpy.int8(1)
    with not_raises(TypeError):
        prob.nobj = numpy.int16(1)
    with not_raises(TypeError):
        prob.nobj = numpy.int32(1)
    with not_raises(TypeError):
        prob.nobj = numpy.int64(1)
    with not_raises(TypeError):
        prob.nobj = int(1)

def test_nobj_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.nobj = None
    with pytest.raises(TypeError):
        prob.nobj = "string"
    with pytest.raises(TypeError):
        prob.nobj = float(1.0)

def test_nobj_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.nobj = int(1)
    with pytest.raises(ValueError):
        prob.nobj = int(0)
    with pytest.raises(ValueError):
        prob.nobj = int(-1)

def test_nobj_fdel(prob):
    del prob.nobj
    with pytest.raises(AttributeError):
        prob.nobj

##############
### obj_wt ###
##############
def test_obj_wt_is_concrete():
    assert_concrete_property(DenseSolution, "obj_wt")

def test_obj_wt_fget(prob):
    assert isinstance(prob.obj_wt, numpy.ndarray)

def test_obj_wt_fset(prob, obj_wt):
    with not_raises(TypeError):
        prob.obj_wt = obj_wt
    with not_raises(TypeError):
        prob.obj_wt = float(1.0)

def test_obj_wt_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.obj_wt = "string"
    with pytest.raises(TypeError):
        prob.obj_wt = None

def test_obj_wt_fset_ValueError(prob, obj_wt):
    with pytest.raises(ValueError):
        prob.obj_wt = numpy.zeros(len(obj_wt)+1, dtype=obj_wt.dtype)

def test_obj_wt_fdel(prob):
    del prob.obj_wt
    with pytest.raises(AttributeError):
        prob.obj_wt

###############
### nineqcv ###
###############
def test_nineqcv_is_concrete():
    assert_concrete_property(DenseSolution, "nineqcv")

def test_nineqcv_fget(prob):
    assert isinstance(prob.nineqcv, Integral)

def test_nineqcv_fset(prob):
    with not_raises(TypeError):
        prob.nineqcv = numpy.int8(1)
    with not_raises(TypeError):
        prob.nineqcv = numpy.int16(1)
    with not_raises(TypeError):
        prob.nineqcv = numpy.int32(1)
    with not_raises(TypeError):
        prob.nineqcv = numpy.int64(1)
    with not_raises(TypeError):
        prob.nineqcv = int(1)

def test_nineqcv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.nineqcv = None
    with pytest.raises(TypeError):
        prob.nineqcv = "string"
    with pytest.raises(TypeError):
        prob.nineqcv = float(1.0)

def test_nineqcv_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.nineqcv = int(1)
    with not_raises(ValueError):
        prob.nineqcv = int(0)
    with pytest.raises(ValueError):
        prob.nineqcv = int(-1)

def test_nineqcv_fdel(prob):
    del prob.nineqcv
    with pytest.raises(AttributeError):
        prob.nineqcv

#################
### ineqcv_wt ###
#################
def test_ineqcv_wt_is_concrete():
    assert_concrete_property(DenseSolution, "ineqcv_wt")

def test_ineqcv_wt_fget(prob):
    assert isinstance(prob.ineqcv_wt, numpy.ndarray)

def test_ineqcv_wt_fset(prob, ineqcv_wt):
    with not_raises(TypeError):
        prob.ineqcv_wt = ineqcv_wt
    with not_raises(TypeError):
        prob.ineqcv_wt = float(1.0)

def test_ineqcv_wt_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ineqcv_wt = "string"
    with pytest.raises(TypeError):
        prob.ineqcv_wt = None

def test_ineqcv_wt_fset_ValueError(prob, ineqcv_wt):
    with pytest.raises(ValueError):
        prob.ineqcv_wt = numpy.zeros(len(ineqcv_wt)+1, dtype=ineqcv_wt.dtype)

def test_ineqcv_wt_fdel(prob):
    del prob.ineqcv_wt
    with pytest.raises(AttributeError):
        prob.ineqcv_wt

#############
### neqcv ###
#############
def test_neqcv_is_concrete():
    assert_concrete_property(DenseSolution, "neqcv")

def test_neqcv_fget(prob):
    assert isinstance(prob.neqcv, Integral)

def test_neqcv_fset(prob):
    with not_raises(TypeError):
        prob.neqcv = numpy.int8(1)
    with not_raises(TypeError):
        prob.neqcv = numpy.int16(1)
    with not_raises(TypeError):
        prob.neqcv = numpy.int32(1)
    with not_raises(TypeError):
        prob.neqcv = numpy.int64(1)
    with not_raises(TypeError):
        prob.neqcv = int(1)

def test_neqcv_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.neqcv = None
    with pytest.raises(TypeError):
        prob.neqcv = "string"
    with pytest.raises(TypeError):
        prob.neqcv = float(1.0)

def test_neqcv_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.neqcv = int(1)
    with not_raises(ValueError):
        prob.neqcv = int(0)
    with pytest.raises(ValueError):
        prob.neqcv = int(-1)

def test_neqcv_fdel(prob):
    del prob.neqcv
    with pytest.raises(AttributeError):
        prob.neqcv

###############
### eqcv_wt ###
###############
def test_eqcv_wt_is_concrete():
    assert_concrete_property(DenseSolution, "eqcv_wt")

def test_eqcv_wt_fget(prob):
    assert isinstance(prob.eqcv_wt, numpy.ndarray)

def test_eqcv_wt_fset(prob, eqcv_wt):
    with not_raises(TypeError):
        prob.eqcv_wt = eqcv_wt
    with not_raises(TypeError):
        prob.eqcv_wt = float(1.0)

def test_eqcv_wt_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.eqcv_wt = "string"
    with pytest.raises(TypeError):
        prob.eqcv_wt = None

def test_eqcv_wt_fset_ValueError(prob, eqcv_wt):
    with pytest.raises(ValueError):
        prob.eqcv_wt = numpy.zeros(len(eqcv_wt)+1, dtype=eqcv_wt.dtype)

def test_eqcv_wt_fdel(prob):
    del prob.eqcv_wt
    with pytest.raises(AttributeError):
        prob.eqcv_wt

###############
### nsoln #####
###############
def test_nsoln_is_concrete():
    assert_concrete_property(DenseSolution, "nsoln")

def test_nsoln_fget(prob):
    assert isinstance(prob.nsoln, Integral)

def test_nsoln_fset(prob, nsoln):
    with not_raises(TypeError):
        prob.nsoln = numpy.int8(nsoln)
    with not_raises(TypeError):
        prob.nsoln = numpy.int16(nsoln)
    with not_raises(TypeError):
        prob.nsoln = numpy.int32(nsoln)
    with not_raises(TypeError):
        prob.nsoln = numpy.int64(nsoln)
    with not_raises(TypeError):
        prob.nsoln = int(nsoln)

def test_nsoln_fset_TypeError(prob, nsoln):
    with pytest.raises(TypeError):
        prob.nsoln = None
    with pytest.raises(TypeError):
        prob.nsoln = "string"
    with pytest.raises(TypeError):
        prob.nsoln = float(nsoln)

def test_nsoln_fset_ValueError(prob, nsoln):
    with not_raises(ValueError):
        prob.nsoln = int(nsoln)
    with not_raises(ValueError):
        prob.nsoln = int(0)
    with pytest.raises(ValueError):
        prob.nsoln = int(-nsoln)

def test_nsoln_fdel(prob):
    del prob.nsoln
    with pytest.raises(AttributeError):
        prob.nsoln

#################
### soln_decn ###
#################
def test_soln_decn_is_concrete():
    assert_concrete_property(DenseSolution, "soln_decn")

def test_soln_decn_fget(prob, nsoln, ndecn):
    assert isinstance(prob.soln_decn, numpy.ndarray)
    assert prob.soln_decn.shape == (nsoln,ndecn)

def test_soln_decn_fset(prob, nsoln, nobj, decn_space_lower, decn_space_upper):
    with not_raises(TypeError):
        prob.soln_decn = numpy.random.randint(decn_space_lower, decn_space_upper, (nsoln,nobj))

def test_soln_decn_fset_TypeError(prob, ndecn):
    with pytest.raises(TypeError):
        prob.soln_decn = None
    with pytest.raises(TypeError):
        prob.soln_decn = "string"
    with pytest.raises(TypeError):
        prob.soln_decn = int(1)
    with pytest.raises(TypeError):
        prob.soln_decn = ndecn * [1]

def test_soln_decn_fset_ValueError(prob, decn_space_lower):
    with pytest.raises(ValueError):
        prob.soln_decn = numpy.array([1])
    with pytest.raises(ValueError):
        prob.soln_decn = decn_space_lower
    with pytest.raises(ValueError):
        prob.soln_decn = numpy.concatenate([decn_space_lower,decn_space_lower])
    with pytest.raises(ValueError):
        prob.soln_decn = numpy.stack([
            decn_space_lower,decn_space_lower,
            decn_space_lower,decn_space_lower
        ])
    with pytest.raises(ValueError):
        prob.soln_decn = numpy.stack([
            numpy.concatenate([decn_space_lower,decn_space_lower]),
            numpy.concatenate([decn_space_lower,decn_space_lower])
        ])

def test_soln_decn_fdel(prob):
    del prob.soln_decn
    with pytest.raises(AttributeError):
        prob.soln_decn

################
### soln_obj ###
################
def test_soln_obj_is_concrete():
    assert_concrete_property(DenseSolution, "soln_obj")

def test_soln_obj_fget(prob, nsoln, nobj):
    assert isinstance(prob.soln_obj, numpy.ndarray)
    assert prob.soln_obj.shape == (nsoln,nobj)

def test_soln_obj_fset(prob, nsoln, nobj):
    with not_raises(TypeError):
        prob.soln_obj = numpy.random.random((nsoln,nobj))

def test_soln_obj_fset_TypeError(prob, nobj):
    with pytest.raises(TypeError):
        prob.soln_obj = None
    with pytest.raises(TypeError):
        prob.soln_obj = "string"
    with pytest.raises(TypeError):
        prob.soln_obj = int(1)
    with pytest.raises(TypeError):
        prob.soln_obj = nobj * [1]

def test_soln_obj_fset_ValueError(prob, nsoln, nobj):
    with pytest.raises(ValueError):
        prob.soln_obj = numpy.random.random(nsoln)
    with pytest.raises(ValueError):
        prob.soln_obj = numpy.random.random((nsoln+1,nobj+1))

def test_soln_obj_fdel(prob):
    del prob.soln_obj
    with pytest.raises(AttributeError):
        prob.soln_obj

###################
### soln_ineqcv ###
###################
def test_soln_ineqcv_is_concrete():
    assert_concrete_property(DenseSolution, "soln_ineqcv")

def test_soln_ineqcv_fget(prob, nsoln, nineqcv):
    assert isinstance(prob.soln_ineqcv, numpy.ndarray)
    assert prob.soln_ineqcv.shape == (nsoln,nineqcv)

def test_soln_ineqcv_fset(prob, nsoln, nineqcv):
    with not_raises(TypeError):
        prob.soln_ineqcv = numpy.random.random((nsoln,nineqcv))

def test_soln_ineqcv_fset_TypeError(prob, nineqcv):
    with pytest.raises(TypeError):
        prob.soln_ineqcv = None
    with pytest.raises(TypeError):
        prob.soln_ineqcv = "string"
    with pytest.raises(TypeError):
        prob.soln_ineqcv = int(1)
    with pytest.raises(TypeError):
        prob.soln_ineqcv = nineqcv * [1]

def test_soln_ineqcv_fset_ValueError(prob, nsoln, nineqcv):
    with pytest.raises(ValueError):
        prob.soln_ineqcv = numpy.random.random(nsoln)
    with pytest.raises(ValueError):
        prob.soln_ineqcv = numpy.random.random((nsoln+1,nineqcv+1))

def test_soln_ineqcv_fdel(prob):
    del prob.soln_ineqcv
    with pytest.raises(AttributeError):
        prob.soln_ineqcv

#################
### soln_eqcv ###
#################
def test_soln_eqcv_is_concrete():
    assert_concrete_property(DenseSolution, "soln_eqcv")

def test_soln_eqcv_fget(prob, nsoln, neqcv):
    assert isinstance(prob.soln_eqcv, numpy.ndarray)
    assert prob.soln_eqcv.shape == (nsoln,neqcv)

def test_soln_eqcv_fset(prob, nsoln, neqcv):
    with not_raises(TypeError):
        prob.soln_eqcv = numpy.random.random((nsoln,neqcv))

def test_soln_eqcv_fset_TypeError(prob, neqcv):
    with pytest.raises(TypeError):
        prob.soln_eqcv = None
    with pytest.raises(TypeError):
        prob.soln_eqcv = "string"
    with pytest.raises(TypeError):
        prob.soln_eqcv = int(1)
    with pytest.raises(TypeError):
        prob.soln_eqcv = neqcv * [1]

def test_soln_eqcv_fset_ValueError(prob, nsoln, neqcv):
    with pytest.raises(ValueError):
        prob.soln_eqcv = numpy.random.random(nsoln)
    with pytest.raises(ValueError):
        prob.soln_eqcv = numpy.random.random((nsoln+1,neqcv+1))

def test_soln_eqcv_fdel(prob):
    del prob.soln_eqcv
    with pytest.raises(AttributeError):
        prob.soln_eqcv

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseSolution, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseSolution_is_concrete():
    assert_concrete_function(check_is_DenseSolution)

def test_check_is_DenseSolution(prob):
    with not_raises(TypeError):
        check_is_DenseSolution(prob, "prob")
    with pytest.raises(TypeError):
        check_is_DenseSolution(None, "prob")
