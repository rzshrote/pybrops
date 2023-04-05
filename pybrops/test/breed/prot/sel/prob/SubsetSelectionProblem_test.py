import pytest

from pybrops.test.assert_python import assert_concrete_function, assert_docstring, not_raises
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_property

from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem, check_is_SubsetSelectionProblem

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def prob():
    yield SubsetSelectionProblem()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_SubsetSelectionProblem_docstring():
    assert_docstring(SubsetSelectionProblem)

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SubsetSelectionProblem, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SubsetSelectionProblem_is_concrete():
    assert_concrete_function(check_is_SubsetSelectionProblem)

def test_check_is_SubsetSelectionProblem(prob):
    with not_raises(TypeError):
        check_is_SubsetSelectionProblem(prob, "prob")
    with pytest.raises(TypeError):
        check_is_SubsetSelectionProblem(None, "prob")
