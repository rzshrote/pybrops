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
def test_encode_trans_is_abstract():
    assert_abstract_property(SelectionProblem, "encode_trans")

def test_encode_trans_kwargs_is_abstract():
    assert_abstract_property(SelectionProblem, "encode_trans_kwargs")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_encodefn_is_abstract(prob):
    assert_abstract_method(prob, "encodefn")

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
