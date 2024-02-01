from numbers import Integral, Number
from typing import Callable, Iterable, Sequence
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.SubsetSolution import SubsetSolution, check_is_SubsetSolution

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
    yield SubsetSolution(
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
def test_DenseSubsetSolution_docstring():
    assert_docstring(SubsetSolution)

################################################################################
########################### Test concrete properties ###########################
################################################################################

##################
### decn_space ###
##################
def test_decn_space_is_concrete():
    assert_concrete_property(SubsetSolution, "decn_space")

def test_decn_space_fget(prob, ndecn):
    assert isinstance(prob.decn_space, numpy.ndarray) or (prob.decn_space is None)
    assert prob.decn_space.shape[0] >= ndecn

def test_decn_space_fset(prob, decn_space_lower, decn_space_upper):
    with not_raises(TypeError,ValueError,Exception):
        prob.decn_space = numpy.concatenate([decn_space_lower,decn_space_upper])

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
        prob.decn_space = numpy.stack([
            decn_space_lower,decn_space_lower,
            decn_space_upper,decn_space_upper
        ])
    with pytest.raises(ValueError):
        prob.decn_space = numpy.stack([
            numpy.concatenate([decn_space_lower,decn_space_lower]),
            numpy.concatenate([decn_space_upper,decn_space_upper])
        ])

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SubsetSolution, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseSubsetSolution_is_concrete():
    assert_concrete_function(check_is_SubsetSolution)

def test_check_is_DenseSubsetSolution(prob):
    with not_raises(TypeError):
        check_is_SubsetSolution(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetSolution(None, "prob")
