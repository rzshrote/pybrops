import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.op.mate.MatingOperator import MatingOperator
from pybropt.breed.op.mate.MatingOperator import is_MatingOperator
from pybropt.breed.op.mate.MatingOperator import check_is_MatingOperator
from pybropt.breed.op.mate.MatingOperator import cond_check_is_MatingOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield MatingOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(MatingOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(MatingOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mate_is_abstract():
    generic_assert_abstract_method(MatingOperator, "mate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_MatingOperator_is_concrete():
    generic_assert_concrete_function(is_MatingOperator)

def test_check_is_MatingOperator_is_concrete():
    generic_assert_concrete_function(check_is_MatingOperator)

def test_cond_check_is_MatingOperator_is_concrete():
    generic_assert_concrete_function(cond_check_is_MatingOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_MatingOperator(operator):
    assert is_MatingOperator(operator)

def test_check_is_MatingOperator(operator):
    with not_raises(TypeError):
        check_is_MatingOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_MatingOperator(None, "operator")
