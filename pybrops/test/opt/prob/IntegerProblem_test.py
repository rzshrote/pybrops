import numpy
import pytest
from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_property
from pybrops.opt.prob.IntegerProblem import IntegerProblem, check_is_IntegerProblem

################################################################################
################################ Test fixtures #################################
################################################################################
class IntegerProblemTestClass(IntegerProblem):
    def __init__(
            self, ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        ):
        """NA"""
        super(IntegerProblemTestClass, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt, **kwargs
        )
    ### method required by PyBrOpS interface ###
    def evalfn(self, x, *args, **kwargs):
        """NA"""
        super(IntegerProblemTestClass, self).evalfn(x, *args, **kwargs)
    ### method required by PyMOO interface ###
    def _evaluate(self, x, out, *args, **kwargs):
        """NA"""
        super(IntegerProblemTestClass, self)._evaluate(x, out, *args, **kwargs)

@pytest.fixture
def ndecn():
    yield 4

@pytest.fixture
def decn_space_lower():
    yield numpy.array([1,2,3,4], dtype=int)

@pytest.fixture
def decn_space_upper():
    yield numpy.array([5,6,7,8], dtype=int)

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
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt
    ):
    yield IntegerProblemTestClass(
        ndecn, decn_space, decn_space_lower, decn_space_upper, 
        nobj, obj_wt, nineqcv, ineqcv_wt, neqcv, eqcv_wt
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_IntegerProblem_docstring():
    assert_docstring(IntegerProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

############################################################
######## properties inherited from pybrops Problem #########
############################################################

##################
### decn_space ###
##################
def test_decn_space_is_concrete():
    assert_concrete_property(IntegerProblem, "decn_space")

def test_decn_space_fget(prob, ndecn):
    assert isinstance(prob.decn_space, numpy.ndarray) or (prob.decn_space is None)
    assert prob.decn_space.shape == (2,ndecn)

def test_decn_space_fset(prob, decn_space_lower, decn_space_upper):
    with not_raises(TypeError,ValueError,Exception):
        prob.decn_space = numpy.stack([decn_space_lower,decn_space_upper])
    with not_raises(TypeError,ValueError,Exception):
        prob.decn_space = None

def test_decn_space_fset_TypeError(prob, ndecn):
    with pytest.raises(TypeError):
        prob.decn_space = object()
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
            numpy.concatenate([decn_space_lower,decn_space_lower]),
            numpy.concatenate([decn_space_upper,decn_space_upper])
        ])
    with pytest.raises(ValueError):
        prob.decn_space = numpy.concatenate([decn_space_lower,decn_space_upper])

def test_decn_space_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.decn_space

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(IntegerProblem, "__init__")

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
def test_check_is_IntegerProblem_is_concrete():
    assert_concrete_function(check_is_IntegerProblem)

def test_check_is_IntegerProblem(prob):
    with not_raises(TypeError):
        check_is_IntegerProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_IntegerProblem(None, "prob")
