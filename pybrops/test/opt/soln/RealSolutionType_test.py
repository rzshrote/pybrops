import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.RealSolutionType import RealSolutionType, check_is_RealSolutionType

################################################################################
################################ Test fixtures #################################
################################################################################

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(RealSolutionType)

################################################################################
############################# Test concrete methods ############################
################################################################################

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
    assert_concrete_function(check_is_RealSolutionType)

def test_check_is_RealSolution():
    with pytest.raises(TypeError):
        check_is_RealSolutionType(None, "prob")
