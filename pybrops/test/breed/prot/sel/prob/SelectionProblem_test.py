from typing import Callable
import numpy
import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem, check_is_SelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
class SelectionProblemTestClass(SelectionProblem):
    def __init__(
            self,
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, obj_trans, obj_trans_kwargs,
            nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs,
            neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs,
            **kwargs
        ):
        """NA"""
        super(SelectionProblemTestClass, self).__init__(
            ndecn, decn_space, decn_space_lower, decn_space_upper, 
            nobj, obj_wt, obj_trans, obj_trans_kwargs,
            nineqcv, ineqcv_wt, ineqcv_trans, ineqcv_trans_kwargs,
            neqcv, eqcv_wt, eqcv_trans, eqcv_trans_kwargs,
            **kwargs
        )
    @property
    def nlatent(self) -> int:
        """nlatent."""
        return 0
    def latentfn(self, x, *args, **kwargs):
        """NA"""
        super(SelectionProblemTestClass, self).latentfn(x, *args, **kwargs)

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
def obj_trans():
    yield None

@pytest.fixture
def obj_trans_kwargs():
    yield None

@pytest.fixture
def nineqcv():
    yield 0

@pytest.fixture
def ineqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def ineqcv_trans():
    yield None

@pytest.fixture
def ineqcv_trans_kwargs():
    yield None

@pytest.fixture
def neqcv():
    yield 0

@pytest.fixture
def eqcv_wt():
    yield numpy.array([], dtype=float)

@pytest.fixture
def eqcv_trans():
    yield None

@pytest.fixture
def eqcv_trans_kwargs():
    yield None

@pytest.fixture
def prob(
    ndecn,
    decn_space,
    decn_space_lower,
    decn_space_upper,
    nobj,
    obj_wt,
    obj_trans,
    obj_trans_kwargs,
    nineqcv,
    ineqcv_wt,
    ineqcv_trans,
    ineqcv_trans_kwargs,
    neqcv,
    eqcv_wt,
    eqcv_trans,
    eqcv_trans_kwargs
):
    yield SelectionProblemTestClass(
        ndecn = ndecn,
        decn_space = decn_space,
        decn_space_lower = decn_space_lower,
        decn_space_upper = decn_space_upper,
        nobj = nobj,
        obj_wt = obj_wt,
        obj_trans = obj_trans,
        obj_trans_kwargs = obj_trans_kwargs,
        nineqcv = nineqcv,
        ineqcv_wt = ineqcv_wt,
        ineqcv_trans = ineqcv_trans,
        ineqcv_trans_kwargs = ineqcv_trans_kwargs,
        neqcv = neqcv,
        eqcv_wt = eqcv_wt,
        eqcv_trans = eqcv_trans,
        eqcv_trans_kwargs = eqcv_trans_kwargs
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SelectionProblem_docstring():
    assert_docstring(SelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

####################
### obj_trans ###
####################
def test_obj_trans_is_concrete():
    assert_concrete_property(SelectionProblem, "obj_trans")

def test_obj_trans_fget(prob):
    assert isinstance(prob.obj_trans, Callable)

def test_obj_trans_fset(prob, obj_trans):
    with not_raises(Exception):
        prob.obj_trans = obj_trans

def test_obj_trans_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.obj_trans = "string"
    with pytest.raises(TypeError):
        prob.obj_trans = []
    with pytest.raises(TypeError):
        prob.obj_trans = {}

def test_obj_trans_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.obj_trans

###########################
### obj_trans_kwargs ###
###########################
def test_obj_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProblem, "obj_trans_kwargs")

def test_obj_trans_kwargs_fget(prob):
    assert isinstance(prob.obj_trans_kwargs, dict)

def test_obj_trans_kwargs_fset(prob):
    with not_raises(Exception):
        prob.obj_trans_kwargs = {"a":1}
    with not_raises(Exception):
        prob.obj_trans_kwargs = None

def test_obj_trans_kwargs_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.obj_trans_kwargs = "string"
    with pytest.raises(TypeError):
        prob.obj_trans_kwargs = []
    with pytest.raises(TypeError):
        prob.obj_trans_kwargs = ()

def test_obj_trans_kwargs_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.obj_trans_kwargs

####################
### ineqcv_trans ###
####################
def test_ineqcv_trans_is_concrete():
    assert_concrete_property(SelectionProblem, "ineqcv_trans")

def test_ineqcv_trans_fget(prob):
    assert isinstance(prob.ineqcv_trans, Callable)

def test_ineqcv_trans_fset(prob, ineqcv_trans):
    with not_raises(Exception):
        prob.ineqcv_trans = ineqcv_trans

def test_ineqcv_trans_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ineqcv_trans = "string"
    with pytest.raises(TypeError):
        prob.ineqcv_trans = []
    with pytest.raises(TypeError):
        prob.ineqcv_trans = {}

def test_ineqcv_trans_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.ineqcv_trans

###########################
### ineqcv_trans_kwargs ###
###########################
def test_ineqcv_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProblem, "ineqcv_trans_kwargs")

def test_ineqcv_trans_kwargs_fget(prob):
    assert isinstance(prob.ineqcv_trans_kwargs, dict)

def test_ineqcv_trans_kwargs_fset(prob):
    with not_raises(Exception):
        prob.ineqcv_trans_kwargs = {"a":1}
    with not_raises(Exception):
        prob.ineqcv_trans_kwargs = None

def test_ineqcv_trans_kwargs_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.ineqcv_trans_kwargs = "string"
    with pytest.raises(TypeError):
        prob.ineqcv_trans_kwargs = []
    with pytest.raises(TypeError):
        prob.ineqcv_trans_kwargs = ()

def test_ineqcv_trans_kwargs_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.ineqcv_trans_kwargs

####################
### eqcv_trans ###
####################
def test_eqcv_trans_is_concrete():
    assert_concrete_property(SelectionProblem, "eqcv_trans")

def test_eqcv_trans_fget(prob):
    assert isinstance(prob.eqcv_trans, Callable)

def test_eqcv_trans_fset(prob, eqcv_trans):
    with not_raises(Exception):
        prob.eqcv_trans = eqcv_trans

def test_eqcv_trans_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.eqcv_trans = "string"
    with pytest.raises(TypeError):
        prob.eqcv_trans = []
    with pytest.raises(TypeError):
        prob.eqcv_trans = {}

def test_eqcv_trans_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.eqcv_trans

###########################
### eqcv_trans_kwargs ###
###########################
def test_eqcv_trans_kwargs_is_concrete():
    assert_concrete_property(SelectionProblem, "eqcv_trans_kwargs")

def test_eqcv_trans_kwargs_fget(prob):
    assert isinstance(prob.eqcv_trans_kwargs, dict)

def test_eqcv_trans_kwargs_fset(prob):
    with not_raises(Exception):
        prob.eqcv_trans_kwargs = {"a":1}
    with not_raises(Exception):
        prob.eqcv_trans_kwargs = None

def test_eqcv_trans_kwargs_fset_TypeError(prob):
    with pytest.raises(TypeError):
        prob.eqcv_trans_kwargs = "string"
    with pytest.raises(TypeError):
        prob.eqcv_trans_kwargs = []
    with pytest.raises(TypeError):
        prob.eqcv_trans_kwargs = ()

def test_eqcv_trans_kwargs_fdel(prob):
    with pytest.raises(AttributeError):
        del prob.eqcv_trans_kwargs

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SelectionProblem, "__init__")

def test_evalfn_is_concrete():
    assert_concrete_method(SelectionProblem, "evalfn")

def test__evaluate_is_concrete():
    assert_concrete_method(SelectionProblem, "_evaluate")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_latentfn_is_abstract(prob):
    assert_abstract_method(prob, "latentfn")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SelectionProblem_is_concrete():
    assert_concrete_function(check_is_SelectionProblem)

def test_check_is_SelectionProblem(prob):
    with not_raises(TypeError):
        check_is_SelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SelectionProblem(None, "prob")
