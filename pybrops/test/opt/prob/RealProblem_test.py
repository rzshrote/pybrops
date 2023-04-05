import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method

from pybrops.opt.prob.RealProblem import RealProblem, check_is_RealProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield RealProblem()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_RealProblem_docstring():
    assert_docstring(RealProblem)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(RealProblem, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_RealProblem_is_concrete():
    assert_concrete_function(check_is_RealProblem)

def test_check_is_RealProblem(prob):
    with not_raises(TypeError):
        check_is_RealProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_RealProblem(None, "prob")
