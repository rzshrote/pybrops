import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.SubsetSolution import SubsetSolution, check_is_SubsetSolution

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SubsetSolution()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(SubsetSolution)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SubsetSolution, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SubsetSolution_is_concrete():
    assert_concrete_function(check_is_SubsetSolution)

def test_check_is_SubsetSolution(prob):
    with not_raises(TypeError):
        check_is_SubsetSolution(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetSolution(None, "prob")
