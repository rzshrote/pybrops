import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.SubsetSelectionProblemType import SubsetSelectionProblemType, check_is_SubsetSelectionProblemType

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SubsetSelectionProblemType()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetSelectionProblemType_docstring():
    assert_docstring(SubsetSelectionProblemType)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SubsetSelectionProblemType, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SubsetSelectionProblemType_is_concrete():
    assert_concrete_function(check_is_SubsetSelectionProblemType)

def test_check_is_SubsetSelectionProblemType(prob):
    with not_raises(TypeError):
        check_is_SubsetSelectionProblemType(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetSelectionProblemType(None, "prob")
