import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.arch.EmigrationOperator import EmigrationOperator
from pybropt.breed.arch.EmigrationOperator import is_EmigrationOperator
from pybropt.breed.arch.EmigrationOperator import check_is_EmigrationOperator
from pybropt.breed.arch.EmigrationOperator import cond_check_is_EmigrationOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield EmigrationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(EmigrationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(EmigrationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(EmigrationOperator, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_EmigrationOperator_is_concrete():
    generic_assert_concrete_function(is_EmigrationOperator)

def test_check_is_EmigrationOperator_is_concrete():
    generic_assert_concrete_function(check_is_EmigrationOperator)

def test_cond_check_is_EmigrationOperator_is_concrete():
    generic_assert_concrete_function(cond_check_is_EmigrationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_EmigrationOperator(arch):
    assert is_EmigrationOperator(arch)

def test_check_is_EmigrationOperator(arch):
    with not_raises(TypeError):
        check_is_EmigrationOperator(arch, "arch")
    with pytest.raises(TypeError):
        check_is_EmigrationOperator(None, "arch")
