import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method

from pybrops.opt.prob.SubsetProblemType import SubsetProblemType, check_is_SubsetProblemType

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SubsetProblemType()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetProblem_docstring():
    assert_docstring(SubsetProblemType)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SubsetProblemType, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SubsetProblem_is_concrete():
    assert_concrete_function(check_is_SubsetProblemType)

def test_check_is_SubsetProblem(prob):
    with not_raises(TypeError):
        check_is_SubsetProblemType(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetProblemType(None, "prob")
