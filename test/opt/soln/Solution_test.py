from numbers import Integral, Number
from typing import Callable, Iterable, Sequence
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.Solution import Solution, check_is_Solution

################################################################################
################################ Test fixtures #################################
################################################################################
class DummySolution(Solution):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, 
            nsoln, soln_decn, soln_obj, soln_ineqcv, soln_eqcv, **kwargs
        ):
        """NA"""
        self.ndecn = ndecn
        self.decn_space = decn_space
        self.decn_space_lower = decn_space_lower
        self.decn_space_upper = decn_space_upper
        self.nobj = nobj
        self.obj_wt = obj_wt
        self.nineqcv = nineqcv
        self.ineqcv_wt = ineqcv_wt
        self.neqcv = neqcv
        self.eqcv_wt = eqcv_wt
        self.nsoln = nsoln
        self.soln_decn = soln_decn
        self.soln_obj = soln_obj
        self.soln_ineqcv = soln_ineqcv
        self.soln_eqcv = soln_eqcv

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
def soln(
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
    yield DummySolution(
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
def test_Solution_docstring():
    assert_docstring(Solution)

################################################################################
########################### Test concrete properties ###########################
################################################################################

#############
### ndecn ###
#############
def test_ndecn_is_concrete():
    assert_concrete_property(Solution, "ndecn")

def test_ndecn_fget(soln):
    assert isinstance(soln.ndecn, Integral)

def test_ndecn_fset(soln):
    with not_raises(TypeError):
        soln.ndecn = numpy.int8(1)
    with not_raises(TypeError):
        soln.ndecn = numpy.int16(1)
    with not_raises(TypeError):
        soln.ndecn = numpy.int32(1)
    with not_raises(TypeError):
        soln.ndecn = numpy.int64(1)
    with not_raises(TypeError):
        soln.ndecn = int(1)

def test_ndecn_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.ndecn = None
    with pytest.raises(TypeError):
        soln.ndecn = "string"
    with pytest.raises(TypeError):
        soln.ndecn = float(1.0)

def test_ndecn_fset_ValueError(soln):
    with pytest.raises(ValueError):
        soln.ndecn = int(0)
    with pytest.raises(ValueError):
        soln.ndecn = int(-1)

##################
### decn_space ###
##################
def test_decn_space_is_concrete():
    assert_concrete_property(Solution, "decn_space")

def test_decn_space_fget(soln, ndecn):
    assert isinstance(soln.decn_space, numpy.ndarray) or (soln.decn_space is None)
    if soln.decn_space is not None:
        assert soln.decn_space.shape == (2,ndecn)

def test_decn_space_fset(soln, decn_space_lower, decn_space_upper):
    with not_raises(TypeError):
        soln.decn_space = numpy.stack([decn_space_lower,decn_space_upper])
    with not_raises(TypeError):
        soln.decn_space = None

def test_decn_space_fset_TypeError(soln, ndecn):
    with pytest.raises(TypeError):
        soln.decn_space = "string"
    with pytest.raises(TypeError):
        soln.decn_space = int(1)
    with pytest.raises(TypeError):
        soln.decn_space = ndecn * [1]

def test_decn_space_fset_ValueError(soln, decn_space_lower, decn_space_upper):
    with pytest.raises(ValueError):
        soln.decn_space = numpy.array([1])
    with pytest.raises(ValueError):
        soln.decn_space = decn_space_lower
    with pytest.raises(ValueError):
        soln.decn_space = numpy.concatenate([decn_space_lower,decn_space_upper])
    with pytest.raises(ValueError):
        soln.decn_space = numpy.stack([
            decn_space_lower,decn_space_lower,
            decn_space_upper,decn_space_upper
        ])
    with pytest.raises(ValueError):
        soln.decn_space = numpy.stack([
            numpy.concatenate([decn_space_lower,decn_space_lower]),
            numpy.concatenate([decn_space_upper,decn_space_upper])
        ])

########################
### decn_space_lower ###
########################
def test_decn_space_lower_is_concrete():
    assert_concrete_property(Solution, "decn_space_lower")

def test_decn_space_lower_fget(soln):
    assert isinstance(soln.decn_space_lower, numpy.ndarray) or (soln.decn_space_lower is None)
    assert len(soln.decn_space_lower) == soln.ndecn

def test_decn_space_lower_fset(soln, decn_space_lower):
    with not_raises(TypeError):
        soln.decn_space_lower = decn_space_lower
    with not_raises(TypeError):
        soln.decn_space_lower = decn_space_lower[0]
    with not_raises(TypeError):
        soln.decn_space_lower = None

def test_decn_space_lower_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.decn_space_lower = "string"

def test_decn_space_lower_fset_ValueError(soln, decn_space_lower):
    with pytest.raises(ValueError):
        soln.decn_space_lower = numpy.concatenate([decn_space_lower, decn_space_lower])

########################
### decn_space_upper ###
########################
def test_decn_space_upper_is_concrete():
    assert_concrete_property(Solution, "decn_space_upper")

def test_decn_space_upper_fget(soln):
    assert isinstance(soln.decn_space_upper, numpy.ndarray) or (soln.decn_space_upper is None)
    assert len(soln.decn_space_upper) == soln.ndecn

def test_decn_space_upper_fset(soln, decn_space_upper):
    with not_raises(TypeError):
        soln.decn_space_upper = decn_space_upper
    with not_raises(TypeError):
        soln.decn_space_upper = decn_space_upper[0]
    with not_raises(TypeError):
        soln.decn_space_upper = None

def test_decn_space_upper_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.decn_space_upper = "string"

def test_decn_space_upper_fset_ValueError(soln, decn_space_upper):
    with pytest.raises(ValueError):
        soln.decn_space_upper = numpy.concatenate([decn_space_upper, decn_space_upper])

############
### nobj ###
############
def test_nobj_is_concrete():
    assert_concrete_property(Solution, "nobj")

def test_nobj_fget(soln):
    assert isinstance(soln.nobj, Integral)

def test_nobj_fset(soln):
    with not_raises(TypeError):
        soln.nobj = numpy.int8(1)
    with not_raises(TypeError):
        soln.nobj = numpy.int16(1)
    with not_raises(TypeError):
        soln.nobj = numpy.int32(1)
    with not_raises(TypeError):
        soln.nobj = numpy.int64(1)
    with not_raises(TypeError):
        soln.nobj = int(1)

def test_nobj_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.nobj = None
    with pytest.raises(TypeError):
        soln.nobj = "string"
    with pytest.raises(TypeError):
        soln.nobj = float(1.0)

def test_nobj_fset_ValueError(soln):
    with not_raises(ValueError):
        soln.nobj = int(1)
    with pytest.raises(ValueError):
        soln.nobj = int(0)
    with pytest.raises(ValueError):
        soln.nobj = int(-1)

##############
### obj_wt ###
##############
def test_obj_wt_is_concrete():
    assert_concrete_property(Solution, "obj_wt")

def test_obj_wt_fget(soln):
    assert isinstance(soln.obj_wt, numpy.ndarray)

def test_obj_wt_fset(soln, obj_wt):
    with not_raises(TypeError):
        soln.obj_wt = obj_wt
    with not_raises(TypeError):
        soln.obj_wt = float(1.0)
    with not_raises(TypeError):
        soln.obj_wt = None

def test_obj_wt_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.obj_wt = object()
    with pytest.raises(TypeError):
        soln.obj_wt = "string"

def test_obj_wt_fset_ValueError(soln, obj_wt):
    with pytest.raises(ValueError):
        soln.obj_wt = numpy.zeros(len(obj_wt)+1, dtype=obj_wt.dtype)

###############
### nineqcv ###
###############
def test_nineqcv_is_concrete():
    assert_concrete_property(Solution, "nineqcv")

def test_nineqcv_fget(soln):
    assert isinstance(soln.nineqcv, Integral)

def test_nineqcv_fset(soln):
    with not_raises(TypeError):
        soln.nineqcv = numpy.int8(1)
    with not_raises(TypeError):
        soln.nineqcv = numpy.int16(1)
    with not_raises(TypeError):
        soln.nineqcv = numpy.int32(1)
    with not_raises(TypeError):
        soln.nineqcv = numpy.int64(1)
    with not_raises(TypeError):
        soln.nineqcv = int(1)
    with not_raises(TypeError):
        soln.nineqcv = None

def test_nineqcv_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.nineqcv = object()
    with pytest.raises(TypeError):
        soln.nineqcv = "string"
    with pytest.raises(TypeError):
        soln.nineqcv = float(1.0)

def test_nineqcv_fset_ValueError(soln):
    with not_raises(ValueError):
        soln.nineqcv = int(1)
    with not_raises(ValueError):
        soln.nineqcv = int(0)
    with pytest.raises(ValueError):
        soln.nineqcv = int(-1)

#################
### ineqcv_wt ###
#################
def test_ineqcv_wt_is_concrete():
    assert_concrete_property(Solution, "ineqcv_wt")

def test_ineqcv_wt_fget(soln):
    assert isinstance(soln.ineqcv_wt, numpy.ndarray)

def test_ineqcv_wt_fset(soln, ineqcv_wt):
    with not_raises(TypeError):
        soln.ineqcv_wt = ineqcv_wt
    with not_raises(TypeError):
        soln.ineqcv_wt = float(1.0)
    with not_raises(TypeError):
        soln.ineqcv_wt = None

def test_ineqcv_wt_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.ineqcv_wt = "string"

def test_ineqcv_wt_fset_ValueError(soln, ineqcv_wt):
    with pytest.raises(ValueError):
        soln.ineqcv_wt = numpy.zeros(len(ineqcv_wt)+1, dtype=ineqcv_wt.dtype)

#############
### neqcv ###
#############
def test_neqcv_is_concrete():
    assert_concrete_property(Solution, "neqcv")

def test_neqcv_fget(soln):
    assert isinstance(soln.neqcv, Integral)

def test_neqcv_fset(soln):
    with not_raises(TypeError):
        soln.neqcv = numpy.int8(1)
    with not_raises(TypeError):
        soln.neqcv = numpy.int16(1)
    with not_raises(TypeError):
        soln.neqcv = numpy.int32(1)
    with not_raises(TypeError):
        soln.neqcv = numpy.int64(1)
    with not_raises(TypeError):
        soln.neqcv = int(1)
    with not_raises(TypeError):
        soln.neqcv = None

def test_neqcv_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.neqcv = "string"
    with pytest.raises(TypeError):
        soln.neqcv = float(1.0)

def test_neqcv_fset_ValueError(soln):
    with not_raises(ValueError):
        soln.neqcv = int(1)
    with not_raises(ValueError):
        soln.neqcv = int(0)
    with pytest.raises(ValueError):
        soln.neqcv = int(-1)

###############
### eqcv_wt ###
###############
def test_eqcv_wt_is_concrete():
    assert_concrete_property(Solution, "eqcv_wt")

def test_eqcv_wt_fget(soln):
    assert isinstance(soln.eqcv_wt, numpy.ndarray)

def test_eqcv_wt_fset(soln, eqcv_wt):
    with not_raises(TypeError):
        soln.eqcv_wt = eqcv_wt
    with not_raises(TypeError):
        soln.eqcv_wt = float(1.0)
    with not_raises(TypeError):
        soln.eqcv_wt = None

def test_eqcv_wt_fset_TypeError(soln):
    with pytest.raises(TypeError):
        soln.eqcv_wt = object()
    with pytest.raises(TypeError):
        soln.eqcv_wt = "string"

def test_eqcv_wt_fset_ValueError(soln, eqcv_wt):
    with pytest.raises(ValueError):
        soln.eqcv_wt = numpy.zeros(len(eqcv_wt)+1, dtype=eqcv_wt.dtype)

###############
### nsoln #####
###############
def test_nsoln_is_concrete():
    assert_concrete_property(Solution, "nsoln")

def test_nsoln_fget(soln):
    assert isinstance(soln.nsoln, Integral)

def test_nsoln_fset(soln, nsoln):
    with not_raises(TypeError):
        soln.nsoln = numpy.int8(nsoln)
    with not_raises(TypeError):
        soln.nsoln = numpy.int16(nsoln)
    with not_raises(TypeError):
        soln.nsoln = numpy.int32(nsoln)
    with not_raises(TypeError):
        soln.nsoln = numpy.int64(nsoln)
    with not_raises(TypeError):
        soln.nsoln = int(nsoln)

def test_nsoln_fset_TypeError(soln, nsoln):
    with pytest.raises(TypeError):
        soln.nsoln = None
    with pytest.raises(TypeError):
        soln.nsoln = "string"
    with pytest.raises(TypeError):
        soln.nsoln = float(nsoln)

def test_nsoln_fset_ValueError(soln, nsoln):
    with not_raises(ValueError):
        soln.nsoln = int(nsoln)
    with not_raises(ValueError):
        soln.nsoln = int(0)
    with pytest.raises(ValueError):
        soln.nsoln = int(-nsoln)

#################
### soln_decn ###
#################
def test_soln_decn_is_concrete():
    assert_concrete_property(Solution, "soln_decn")

def test_soln_decn_fget(soln, nsoln, ndecn):
    assert isinstance(soln.soln_decn, numpy.ndarray)
    assert soln.soln_decn.shape == (nsoln,ndecn)

def test_soln_decn_fset(soln, nsoln, nobj, decn_space_lower, decn_space_upper):
    with not_raises(TypeError):
        soln.soln_decn = numpy.random.randint(decn_space_lower, decn_space_upper, (nsoln,nobj))

def test_soln_decn_fset_TypeError(soln, ndecn):
    with pytest.raises(TypeError):
        soln.soln_decn = None
    with pytest.raises(TypeError):
        soln.soln_decn = "string"
    with pytest.raises(TypeError):
        soln.soln_decn = int(1)
    with pytest.raises(TypeError):
        soln.soln_decn = ndecn * [1]

def test_soln_decn_fset_ValueError(soln, decn_space_lower):
    with pytest.raises(ValueError):
        soln.soln_decn = numpy.array([1])
    with pytest.raises(ValueError):
        soln.soln_decn = decn_space_lower
    with pytest.raises(ValueError):
        soln.soln_decn = numpy.concatenate([decn_space_lower,decn_space_lower])
    with pytest.raises(ValueError):
        soln.soln_decn = numpy.stack([
            decn_space_lower,decn_space_lower,
            decn_space_lower,decn_space_lower
        ])
    with pytest.raises(ValueError):
        soln.soln_decn = numpy.stack([
            numpy.concatenate([decn_space_lower,decn_space_lower]),
            numpy.concatenate([decn_space_lower,decn_space_lower])
        ])

################
### soln_obj ###
################
def test_soln_obj_is_concrete():
    assert_concrete_property(Solution, "soln_obj")

def test_soln_obj_fget(soln, nsoln, nobj):
    assert isinstance(soln.soln_obj, numpy.ndarray)
    assert soln.soln_obj.shape == (nsoln,nobj)

def test_soln_obj_fset(soln, nsoln, nobj):
    with not_raises(TypeError):
        soln.soln_obj = numpy.random.random((nsoln,nobj))

def test_soln_obj_fset_TypeError(soln, nobj):
    with pytest.raises(TypeError):
        soln.soln_obj = None
    with pytest.raises(TypeError):
        soln.soln_obj = "string"
    with pytest.raises(TypeError):
        soln.soln_obj = int(1)
    with pytest.raises(TypeError):
        soln.soln_obj = nobj * [1]

def test_soln_obj_fset_ValueError(soln, nsoln, nobj):
    with pytest.raises(ValueError):
        soln.soln_obj = numpy.random.random(nsoln)
    with pytest.raises(ValueError):
        soln.soln_obj = numpy.random.random((nsoln+1,nobj+1))

###################
### soln_ineqcv ###
###################
def test_soln_ineqcv_is_concrete():
    assert_concrete_property(Solution, "soln_ineqcv")

def test_soln_ineqcv_fget(soln, nsoln, nineqcv):
    assert isinstance(soln.soln_ineqcv, numpy.ndarray)
    assert soln.soln_ineqcv.shape == (nsoln,nineqcv)

def test_soln_ineqcv_fset(soln, nsoln, nineqcv):
    with not_raises(TypeError):
        soln.soln_ineqcv = numpy.random.random((nsoln,nineqcv))
    with not_raises(TypeError):
        soln.soln_ineqcv = None

def test_soln_ineqcv_fset_TypeError(soln, nineqcv):
    with pytest.raises(TypeError):
        soln.soln_ineqcv = object()
    with pytest.raises(TypeError):
        soln.soln_ineqcv = "string"
    with pytest.raises(TypeError):
        soln.soln_ineqcv = int(1)
    with pytest.raises(TypeError):
        soln.soln_ineqcv = nineqcv * [1]

def test_soln_ineqcv_fset_ValueError(soln, nsoln, nineqcv):
    with pytest.raises(ValueError):
        soln.soln_ineqcv = numpy.random.random(nsoln)
    with pytest.raises(ValueError):
        soln.soln_ineqcv = numpy.random.random((nsoln+1,nineqcv+1))

#################
### soln_eqcv ###
#################
def test_soln_eqcv_is_concrete():
    assert_concrete_property(Solution, "soln_eqcv")

def test_soln_eqcv_fget(soln, nsoln, neqcv):
    assert isinstance(soln.soln_eqcv, numpy.ndarray)
    assert soln.soln_eqcv.shape == (nsoln,neqcv)

def test_soln_eqcv_fset(soln, nsoln, neqcv):
    with not_raises(TypeError):
        soln.soln_eqcv = numpy.random.random((nsoln,neqcv))
    with not_raises(TypeError):
        soln.soln_eqcv = None

def test_soln_eqcv_fset_TypeError(soln, neqcv):
    with pytest.raises(TypeError):
        soln.soln_eqcv = object()
    with pytest.raises(TypeError):
        soln.soln_eqcv = "string"
    with pytest.raises(TypeError):
        soln.soln_eqcv = int(1)
    with pytest.raises(TypeError):
        soln.soln_eqcv = neqcv * [1]

def test_soln_eqcv_fset_ValueError(soln, nsoln, neqcv):
    with pytest.raises(ValueError):
        soln.soln_eqcv = numpy.random.random(nsoln)
    with pytest.raises(ValueError):
        soln.soln_eqcv = numpy.random.random((nsoln+1,neqcv+1))

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(Solution, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_Solution_is_concrete():
    assert_concrete_function(check_is_Solution)

def test_check_is_Solution(soln):
    with not_raises(TypeError):
        check_is_Solution(soln, "soln")
    with pytest.raises(TypeError):
        check_is_Solution(None, "soln")
