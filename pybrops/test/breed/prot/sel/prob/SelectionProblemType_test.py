import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.SelectionProblemType import SelectionProblemType, check_is_SelectionProblemType

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SelectionProblemType()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SelectionProblem_docstring():
    assert_docstring(SelectionProblemType)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SelectionProblemType, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nlatent_is_abstract():
    assert_abstract_property(SelectionProblemType, "nlatent")

def test_obj_trans_is_abstract():
    assert_abstract_property(SelectionProblemType, "obj_trans")

def test_obj_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblemType, "obj_trans_kwargs")

def test_ineqcv_trans_is_abstract():
    assert_abstract_property(SelectionProblemType, "ineqcv_trans")

def test_ineqcv_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblemType, "ineqcv_trans_kwargs")

def test_eqcv_trans_is_abstract():
    assert_abstract_property(SelectionProblemType, "eqcv_trans")

def test_eqcv_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblemType, "eqcv_trans_kwargs")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_latentfn_is_abstract(prob):
    assert_abstract_method(prob, "latentfn")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SelectionProblem_is_concrete():
    assert_concrete_function(check_is_SelectionProblemType)

def test_check_is_SelectionProblem(prob):
    with not_raises(TypeError):
        check_is_SelectionProblemType(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SelectionProblemType(None, "prob")
