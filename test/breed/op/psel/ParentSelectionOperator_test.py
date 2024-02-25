import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator
from pybrops.breed.op.psel.ParentSelectionOperator import check_is_ParentSelectionOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield DummyParentSelectionOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(ParentSelectionOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(ParentSelectionOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_pselect_is_abstract():
    assert_method_isabstract(ParentSelectionOperator, "pselect")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_ParentSelectionOperator_is_concrete():
    assert_function_isconcrete(check_is_ParentSelectionOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_ParentSelectionOperator(operator):
    with not_raises(TypeError):
        check_is_ParentSelectionOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_ParentSelectionOperator(None, "operator")
