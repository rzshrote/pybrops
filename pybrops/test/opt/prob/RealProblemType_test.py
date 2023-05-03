import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method

from pybrops.opt.prob.RealProblemType import RealProblemType, check_is_RealProblemType

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield RealProblemType()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_RealProblem_docstring():
    assert_docstring(RealProblemType)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(RealProblemType, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_RealProblem_is_concrete():
    assert_concrete_function(check_is_RealProblemType)

def test_check_is_RealProblem(prob):
    with not_raises(TypeError):
        check_is_RealProblemType(prob, "prob")
    with pytest.raises(TypeError):
        check_is_RealProblemType(None, "prob")
