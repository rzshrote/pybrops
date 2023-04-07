import pytest

from pybrops.test.assert_python import assert_abstract_property_fget, assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem, check_is_SelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SelectionProblem()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SelectionProblem_docstring():
    assert_docstring(SelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SelectionProblem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nlatent_is_abstract():
    assert_abstract_property_fget(SelectionProblem, "nlatent")

def test_obj_trans_is_abstract():
    assert_abstract_property(SelectionProblem, "obj_trans")

def test_obj_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblem, "obj_trans_kwargs")

def test_ineqcv_trans_is_abstract():
    assert_abstract_property(SelectionProblem, "ineqcv_trans")

def test_ineqcv_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblem, "ineqcv_trans_kwargs")

def test_eqcv_trans_is_abstract():
    assert_abstract_property(SelectionProblem, "eqcv_trans")

def test_eqcv_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblem, "eqcv_trans_kwargs")

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
