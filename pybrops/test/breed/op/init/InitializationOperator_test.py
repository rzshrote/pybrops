import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.init.InitializationOperator import is_InitializationOperator
from pybrops.breed.op.init.InitializationOperator import check_is_InitializationOperator
from pybrops.breed.op.init.InitializationOperator import cond_check_is_InitializationOperator

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
    generic_assert_docstring(InitializationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(InitializationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_initialize_is_abstract():
    generic_assert_abstract_method(InitializationOperator, "initialize")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_InitializationOperator_is_concrete():
    generic_assert_concrete_function(is_InitializationOperator)

def test_check_is_InitializationOperator_is_concrete():
    generic_assert_concrete_function(check_is_InitializationOperator)

def test_cond_check_is_InitializationOperator_is_concrete():
    generic_assert_concrete_function(cond_check_is_InitializationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_InitializationOperator(operator):
    assert is_InitializationOperator(operator)

def test_check_is_InitializationOperator(operator):
    with not_raises(TypeError):
        check_is_InitializationOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_InitializationOperator(None, "operator")
