from numbers import Integral, Number
from typing import Callable, Iterable, Sequence
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.prob.DenseProblem import DenseProblem, check_is_DenseProblem

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
    eqcv_wt
):
    yield DenseProblem(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj,
        obj_wt = obj_wt,
        nineqcv = nineqcv,
        ineqcv_wt = ineqcv_wt,
        neqcv = neqcv,
        eqcv_wt = eqcv_wt
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_DenseProblem_docstring():
    assert_docstring(DenseProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

############################################################
######### properties inherited from pymoo Problem ##########
############################################################

#############
### n_var ###
#############
def test_n_var_is_concrete():
    assert_concrete_property(DenseProblem, "n_var")

def test_n_var_fget(prob):
    assert isinstance(prob.n_var, Integral)

def test_n_var_fset(prob):
    with not_raises(TypeError):
        prob.n_var = numpy.int8(1)
    with not_raises(TypeError):
        prob.n_var = numpy.int16(1)
    with not_raises(TypeError):
        prob.n_var = numpy.int32(1)
    with not_raises(TypeError):
        prob.n_var = numpy.int64(1)
    with not_raises(TypeError):
        prob.n_var = int(1)

def test_n_var_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.n_var = None
    with pytest.raises(TypeError):
        prob.n_var = "string"
    with pytest.raises(TypeError):
        prob.n_var = float(1.0)

def test_n_var_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.n_var = int(1)
    with not_raises(ValueError):
        prob.n_var = int(0)
    with not_raises(ValueError):
        prob.n_var = int(-1)

def test_n_var_fdel(prob):
    del prob.n_var
    with pytest.raises(AttributeError):
        prob.n_var

#############
### n_obj ###
#############
def test_n_obj_is_concrete():
    assert_concrete_property(DenseProblem, "n_obj")

def test_n_obj_fget(prob):
    assert isinstance(prob.n_obj, Integral)

def test_n_obj_fset(prob):
    with not_raises(TypeError):
        prob.n_obj = numpy.int8(1)
    with not_raises(TypeError):
        prob.n_obj = numpy.int16(1)
    with not_raises(TypeError):
        prob.n_obj = numpy.int32(1)
    with not_raises(TypeError):
        prob.n_obj = numpy.int64(1)
    with not_raises(TypeError):
        prob.n_obj = int(1)

def test_n_obj_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.n_obj = None
    with pytest.raises(TypeError):
        prob.n_obj = "string"
    with pytest.raises(TypeError):
        prob.n_obj = float(1.0)

def test_n_obj_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.n_obj = int(1)
    with pytest.raises(ValueError):
        prob.n_obj = int(0)
    with pytest.raises(ValueError):
        prob.n_obj = int(-1)

def test_n_obj_fdel(prob):
    del prob.n_obj
    with pytest.raises(AttributeError):
        prob.n_obj

####################
### n_ieq_constr ###
####################
def test_n_ieq_constr_is_concrete():
    assert_concrete_property(DenseProblem, "n_ieq_constr")

def test_n_ieq_constr_fget(prob):
    assert isinstance(prob.n_ieq_constr, Integral)

def test_n_ieq_constr_fset(prob):
    with not_raises(TypeError):
        prob.n_ieq_constr = numpy.int8(1)
    with not_raises(TypeError):
        prob.n_ieq_constr = numpy.int16(1)
    with not_raises(TypeError):
        prob.n_ieq_constr = numpy.int32(1)
    with not_raises(TypeError):
        prob.n_ieq_constr = numpy.int64(1)
    with not_raises(TypeError):
        prob.n_ieq_constr = int(1)

def test_n_ieq_constr_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.n_ieq_constr = None
    with pytest.raises(TypeError):
        prob.n_ieq_constr = "string"
    with pytest.raises(TypeError):
        prob.n_ieq_constr = float(1.0)

def test_n_ieq_constr_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.n_ieq_constr = int(1)
    with not_raises(ValueError):
        prob.n_ieq_constr = int(0)
    with pytest.raises(ValueError):
        prob.n_ieq_constr = int(-1)

def test_n_ieq_constr_fdel(prob):
    del prob.n_ieq_constr
    with pytest.raises(AttributeError):
        prob.n_ieq_constr

###################
### n_eq_constr ###
###################
def test_n_eq_constr_is_concrete():
    assert_concrete_property(DenseProblem, "n_eq_constr")

def test_n_eq_constr_fget(prob):
    assert isinstance(prob.n_eq_constr, Integral)

def test_n_eq_constr_fset(prob):
    with not_raises(TypeError):
        prob.n_eq_constr = numpy.int8(1)
    with not_raises(TypeError):
        prob.n_eq_constr = numpy.int16(1)
    with not_raises(TypeError):
        prob.n_eq_constr = numpy.int32(1)
    with not_raises(TypeError):
        prob.n_eq_constr = numpy.int64(1)
    with not_raises(TypeError):
        prob.n_eq_constr = int(1)

def test_n_eq_constr_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.n_eq_constr = None
    with pytest.raises(TypeError):
        prob.n_eq_constr = "string"
    with pytest.raises(TypeError):
        prob.n_eq_constr = float(1.0)

def test_n_eq_constr_fset_ValueError(prob):
    with not_raises(ValueError):
        prob.n_eq_constr = int(1)
    with not_raises(ValueError):
        prob.n_eq_constr = int(0)
    with pytest.raises(ValueError):
        prob.n_eq_constr = int(-1)

def test_n_eq_constr_fdel(prob):
    del prob.n_eq_constr
    with pytest.raises(AttributeError):
        prob.n_eq_constr

##########
### xl ###
##########
def test_xl_is_concrete():
    assert_concrete_property(DenseProblem, "xl")

def test_xl_fget(prob):
    assert isinstance(prob.xl, numpy.ndarray) or (prob.xl is None)
    assert len(prob.xl) == prob.n_var

def test_xl_fset(prob, decn_space_lower):
    with not_raises(TypeError):
        prob.xl = decn_space_lower
    with not_raises(TypeError):
        prob.xl = decn_space_lower[0]
    with not_raises(TypeError):
        prob.xl = None

def test_xl_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.xl = "string"

def test_xl_fset_ValueError(prob, decn_space_lower):
    with pytest.raises(ValueError):
        prob.xl = numpy.concatenate([decn_space_lower, decn_space_lower])

def test_xl_fdel(prob):
    del prob.xl
    with pytest.raises(AttributeError):
        prob.xl

##########
### xu ###
##########
def test_xu_is_concrete():
    assert_concrete_property(DenseProblem, "xu")

def test_xu_fget(prob):
    assert isinstance(prob.xu, numpy.ndarray) or (prob.xu is None)
    assert len(prob.xu) == prob.n_var

def test_xu_fset(prob, decn_space_upper):
    with not_raises(TypeError):
        prob.xu = decn_space_upper
    with not_raises(TypeError):
        prob.xu = decn_space_upper[0]
    with not_raises(TypeError):
        prob.xu = None

def test_xu_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.xu = "string"

def test_xu_fset_ValueError(prob, decn_space_upper):
    with pytest.raises(ValueError):
        prob.xu = numpy.concatenate([decn_space_upper, decn_space_upper])

def test_xu_fdel(prob):
    del prob.xu
    with pytest.raises(AttributeError):
        prob.xu

#############
### vtype ###
#############
def test_vtype_is_concrete():
    assert_concrete_property(DenseProblem, "vtype")

def test_vtype_fget(prob):
    assert isinstance(prob.vtype, type) or (prob.vtype is None)

############
### vars ###
############
def test_vars_is_concrete():
    assert_concrete_property(DenseProblem, "vars")

def test_vars_fget(prob):
    assert isinstance(prob.vars, Sequence) or (prob.vars is None)

def test_vars_fset(prob):
    with not_raises(TypeError):
        prob.vars = numpy.array([4,4,4,4])
    with not_raises(TypeError):
        prob.vars = [4,4,4,4]
    with not_raises(TypeError):
        prob.vars = None

def test_vars_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.vars = int(4)

def test_vars_fdel(prob):
    del prob.vars
    with pytest.raises(AttributeError):
        prob.vars

###################
### elementwise ###
###################
def test_elementwise_is_concrete():
    assert_concrete_property(DenseProblem, "elementwise")

def test_elementwise_fget(prob):
    assert isinstance(prob.elementwise, bool)

def test_elementwise_fset(prob):
    with not_raises(TypeError):
        prob.elementwise = False

def test_elementwise_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.elementwise = None

def test_elementwise_fdel(prob):
    del prob.elementwise
    with pytest.raises(AttributeError):
        prob.elementwise

########################
### elementwise_func ###
########################
def test_elementwise_func_is_concrete():
    assert_concrete_property(DenseProblem, "elementwise_func")

def test_elementwise_func_fget(prob):
    assert isinstance(prob.elementwise_func, type)

def test_elementwise_func_fset(prob):
    with not_raises(TypeError):
        prob.elementwise_func = type

def test_elementwise_func_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.elementwise_func = None

def test_elementwise_func_fdel(prob):
    del prob.elementwise_func
    with pytest.raises(AttributeError):
        prob.elementwise_func

##########################
### elementwise_runner ###
##########################
def test_elementwise_runner_is_concrete():
    assert_concrete_property(DenseProblem, "elementwise_runner")

def test_elementwise_runner_fget(prob):
    assert isinstance(prob.elementwise_runner, Callable)

def test_elementwise_runner_fset(prob):
    def fn(a,b):
        pass
    with not_raises(TypeError):
        prob.elementwise_runner = fn

def test_elementwise_runner_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.elementwise_runner = None

#############################
### replace_nan_values_by ###
#############################
def test_replace_nan_values_by_is_concrete():
    assert_concrete_property(DenseProblem, "replace_nan_values_by")

def test_replace_nan_values_by_fget(prob):
    assert isinstance(prob.replace_nan_values_by, Number) or (prob.replace_nan_values_by is None)

def test_replace_nan_values_by_fset(prob):
    with not_raises(TypeError):
        prob.replace_nan_values_by = 1.0
    with not_raises(TypeError):
        prob.replace_nan_values_by = None

def test_replace_nan_values_by_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.replace_nan_values_by = "string"

def test_replace_nan_values_by_fdel(prob):
    del prob.replace_nan_values_by
    with pytest.raises(AttributeError):
        prob.replace_nan_values_by

##################################
### exclude_from_serialization ###
##################################
def test_exclude_from_serialization_is_concrete():
    assert_concrete_property(DenseProblem, "exclude_from_serialization")

def test_exclude_from_serialization_fget(prob):
    assert isinstance(prob.exclude_from_serialization, Iterable) or (prob.exclude_from_serialization is None)

def test_exclude_from_serialization_fset(prob):
    with not_raises(TypeError):
        prob.exclude_from_serialization = ["a","b","c"]
    with not_raises(TypeError):
        prob.exclude_from_serialization = None

def test_exclude_from_serialization_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.exclude_from_serialization = int(1)

def test_exclude_from_serialization_fdel(prob):
    del prob.exclude_from_serialization
    with pytest.raises(AttributeError):
        prob.exclude_from_serialization

################
### callback ###
################
def test_callback_is_concrete():
    assert_concrete_property(DenseProblem, "callback")

def test_callback_fget(prob):
    assert isinstance(prob.callback, Callable) or (prob.callback is None)

def test_callback_fset(prob):
    def fn(a,b):
        pass
    with not_raises(TypeError):
        prob.callback = fn
    with not_raises(TypeError):
        prob.callback = None

def test_callback_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.callback = int(1)

def test_callback_fdel(prob):
    del prob.callback
    with pytest.raises(AttributeError):
        prob.callback

##############
### strict ###
##############
def test_strict_is_concrete():
    assert_concrete_property(DenseProblem, "strict")

def test_strict_fget(prob):
    assert isinstance(prob.strict, bool)

def test_strict_fset(prob):
    with not_raises(TypeError):
        prob.strict = False

def test_strict_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.strict = "string"
    with pytest.raises(TypeError):
        prob.strict = None

def test_strict_fdel(prob):
    del prob.strict
    with pytest.raises(AttributeError):
        prob.strict

############
### data ###
############
def test_data_is_concrete():
    assert_concrete_property(DenseProblem, "data")

def test_data_fget(prob):
    assert isinstance(prob.data, dict)

def test_data_fset(prob):
    with not_raises(TypeError):
        prob.data = {}

def test_data_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.data = int(1)

def test_data_fdel(prob):
    del prob.data
    with pytest.raises(AttributeError):
        prob.data

############################################################
######## properties inherited from pybrops Problem #########
############################################################

#############
### ndecn ###
#############
def test_ndecn_is_concrete():
    assert_concrete_property(DenseProblem, "ndecn")

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
    with not_raises(ValueError):
        prob.ndecn = int(1)
    with not_raises(ValueError):
        prob.ndecn = int(0)
    with not_raises(ValueError):
        prob.ndecn = int(-1)

def test_ndecn_fdel(prob):
    del prob.ndecn
    with pytest.raises(AttributeError):
        prob.ndecn

##################
### decn_space ###
##################
def test_decn_space_is_concrete():
    assert_concrete_property(DenseProblem, "decn_space")

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
    assert_concrete_property(DenseProblem, "decn_space_lower")

def test_decn_space_lower_fget(prob):
    assert isinstance(prob.decn_space_lower, numpy.ndarray) or (prob.decn_space_lower is None)
    assert len(prob.decn_space_lower) == prob.n_var

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
    assert_concrete_property(DenseProblem, "decn_space_upper")

def test_decn_space_upper_fget(prob):
    assert isinstance(prob.decn_space_upper, numpy.ndarray) or (prob.decn_space_upper is None)
    assert len(prob.decn_space_upper) == prob.n_var

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
    assert_concrete_property(DenseProblem, "nobj")

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
    assert_concrete_property(DenseProblem, "obj_wt")

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
    assert_concrete_property(DenseProblem, "nineqcv")

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
    assert_concrete_property(DenseProblem, "ineqcv_wt")

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
    assert_concrete_property(DenseProblem, "neqcv")

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
    assert_concrete_property(DenseProblem, "eqcv_wt")

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

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseProblem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_evalfn_is_abstract(prob):
    assert_abstract_method(prob, "evalfn")

def test__evaluate_is_abstract(prob):
    assert_abstract_method(prob, "_evaluate")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseProblem_is_concrete():
    assert_concrete_function(check_is_DenseProblem)

def test_check_is_DenseProblem(prob):
    with not_raises(TypeError):
        check_is_DenseProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_DenseProblem(None, "prob")
