from typing import Callable
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.DenseSelectionProblem import DenseSelectionProblem, check_is_DenseSelectionProblem

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
def encode_trans():
    def fn(latent: numpy.ndarray, **kwargs: dict):
        obj = latent
        ineqcv = numpy.array([], dtype=latent.dtype)
        eqcv = numpy.array([], dtype=latent.dtype)
        return (obj, ineqcv, eqcv)
    yield fn

@pytest.fixture
def encode_trans_kwargs():
    yield {}

@pytest.fixture
def prob(
    ndecn,
    decn_space,
    decn_space_lower,
    decn_space_upper,
    encode_trans,
    encode_trans_kwargs,
    nobj,
    obj_wt,
    nineqcv,
    ineqcv_wt,
    neqcv,
    eqcv_wt
):
    yield DenseSelectionProblem(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        encode_trans = encode_trans,
        encode_trans_kwargs = encode_trans_kwargs,
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
def test_DenseSelectionProblem_docstring():
    assert_docstring(DenseSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

####################
### encode_trans ###
####################
def test_encode_trans_is_concrete():
    assert_concrete_property(DenseSelectionProblem, "encode_trans")

def test_encode_trans_fget(prob):
    assert isinstance(prob.encode_trans, Callable)

def test_encode_trans_fset(prob, encode_trans):
    with not_raises(Exception):
        prob.encode_trans = encode_trans

def test_encode_trans_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.encode_trans = "string"
    with pytest.raises(TypeError):
        prob.encode_trans = []
    with pytest.raises(TypeError):
        prob.encode_trans = {}

def test_encode_trans_fdel(prob):
    del prob.encode_trans
    with pytest.raises(AttributeError):
        prob.encode_trans

###########################
### encode_trans_kwargs ###
###########################
def test_encode_trans_kwargs_is_concrete():
    assert_concrete_property(DenseSelectionProblem, "encode_trans_kwargs")

def test_encode_trans_kwargs_fget(prob):
    assert isinstance(prob.encode_trans_kwargs, dict)

def test_encode_trans_kwargs_fset(prob):
    with not_raises(Exception):
        prob.encode_trans_kwargs = {"a":1}
    with not_raises(Exception):
        prob.encode_trans_kwargs = None

def test_encode_trans_kwargs_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.encode_trans_kwargs = "string"
    with pytest.raises(TypeError):
        prob.encode_trans_kwargs = []
    with pytest.raises(TypeError):
        prob.encode_trans_kwargs = ()

def test_encode_trans_kwargs_fdel(prob):
    del prob.encode_trans_kwargs
    with pytest.raises(AttributeError):
        prob.encode_trans_kwargs

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseSelectionProblem, "__init__")

def test_evalfn_is_concrete():
    assert_concrete_method(DenseSelectionProblem, "evalfn")

def test__evaluate_is_concrete():
    assert_concrete_method(DenseSelectionProblem, "_evaluate")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_encodefn_is_abstract(prob):
    assert_abstract_method(prob, "encodefn")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseSelectionProblem_is_concrete():
    assert_concrete_function(check_is_DenseSelectionProblem)

def test_check_is_DenseSelectionProblem(prob):
    with not_raises(TypeError):
        check_is_DenseSelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_DenseSelectionProblem(None, "prob")
