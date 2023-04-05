from typing import Callable
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.DenseSubsetSelectionProblem import DenseSubsetSelectionProblem, check_is_DenseSubsetSelectionProblem

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
    yield numpy.concatenate([decn_space_lower, decn_space_upper])

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
    yield DenseSubsetSelectionProblem(
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
def test_DenseSubsetSelectionProblem_docstring():
    assert_docstring(DenseSubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseSubsetSelectionProblem, "__init__")

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
def test_check_is_DenseSubsetSelectionProblem_is_concrete():
    assert_concrete_function(check_is_DenseSubsetSelectionProblem)

def test_check_is_DenseSubsetSelectionProblem(prob):
    with not_raises(TypeError):
        check_is_DenseSubsetSelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_DenseSubsetSelectionProblem(None, "prob")
