import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.init.InitializationOperator import check_is_InitializationOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield DummyInitializationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(InitializationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_initialize_is_abstract():
    assert_method_isabstract(InitializationOperator, "initialize")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_InitializationOperator_is_concrete():
    assert_function_isconcrete(check_is_InitializationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_InitializationOperator(operator):
    with not_raises(TypeError):
        check_is_InitializationOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_InitializationOperator(None, "operator")
