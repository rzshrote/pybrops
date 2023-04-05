import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.RealSolution import RealSolution, check_is_RealSolution

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield RealSolution()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(RealSolution)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(RealSolution, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_RealSolution_is_concrete():
    assert_concrete_function(check_is_RealSolution)

def test_check_is_RealSolution(prob):
    with not_raises(TypeError):
        check_is_RealSolution(prob, "prob")
    with pytest.raises(TypeError):
        check_is_RealSolution(None, "prob")
