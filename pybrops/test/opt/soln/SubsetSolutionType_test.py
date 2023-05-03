import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.opt.soln.SubsetSolutionType import SubsetSolutionType, check_is_SubsetSolutionType

################################################################################
################################ Test fixtures #################################
################################################################################

################################################################################
############################## Test class docstring ############################
################################################################################
def test_Problem_docstring():
    assert_docstring(SubsetSolutionType)

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
def test_check_is_SubsetSolution_is_concrete():
    assert_concrete_function(check_is_SubsetSolutionType)

def test_check_is_SubsetSolution():
    with pytest.raises(TypeError):
        check_is_SubsetSolutionType(None, "prob")
