import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method

from pybrops.opt.prob.SubsetProblem import SubsetProblem, check_is_SubsetProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SubsetProblem()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetProblem_docstring():
    assert_docstring(SubsetProblem)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SubsetProblem, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SubsetProblem_is_concrete():
    assert_concrete_function(check_is_SubsetProblem)

def test_check_is_SubsetProblem(prob):
    with not_raises(TypeError):
        check_is_SubsetProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetProblem(None, "prob")
