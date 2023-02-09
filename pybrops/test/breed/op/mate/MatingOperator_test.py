import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.mate.MatingOperator import is_MatingOperator
from pybrops.breed.op.mate.MatingOperator import check_is_MatingOperator

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
