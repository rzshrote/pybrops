import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.init.InitializationOperator import check_is_InitializationOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield InitializationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(InitializationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(InitializationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_initialize_is_abstract():
    assert_abstract_method(InitializationOperator, "initialize")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_InitializationOperator_is_concrete():
    assert_concrete_function(check_is_InitializationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_InitializationOperator(operator):
    with not_raises(TypeError):
        check_is_InitializationOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_InitializationOperator(None, "operator")
