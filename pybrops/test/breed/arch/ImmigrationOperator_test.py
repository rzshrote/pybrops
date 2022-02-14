import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.arch.ImmigrationOperator import ImmigrationOperator
from pybrops.breed.arch.ImmigrationOperator import is_ImmigrationOperator
from pybrops.breed.arch.ImmigrationOperator import check_is_ImmigrationOperator
from pybrops.breed.arch.ImmigrationOperator import cond_check_is_ImmigrationOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield ImmigrationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(ImmigrationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(ImmigrationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(ImmigrationOperator, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_ImmigrationOperator_is_concrete():
    generic_assert_concrete_function(is_ImmigrationOperator)

def test_check_is_ImmigrationOperator_is_concrete():
    generic_assert_concrete_function(check_is_ImmigrationOperator)

def test_cond_check_is_ImmigrationOperator_is_concrete():
    generic_assert_concrete_function(cond_check_is_ImmigrationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_ImmigrationOperator(arch):
    assert is_ImmigrationOperator(arch)

def test_check_is_ImmigrationOperator(arch):
    with not_raises(TypeError):
        check_is_ImmigrationOperator(arch, "arch")
    with pytest.raises(TypeError):
        check_is_ImmigrationOperator(None, "arch")
