import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.op.ssel import SurvivorSelectionOperator
from pybropt.breed.op.ssel import is_SurvivorSelectionOperator
from pybropt.breed.op.ssel import check_is_SurvivorSelectionOperator
from pybropt.breed.op.ssel import cond_check_is_SurvivorSelectionOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield SurvivorSelectionOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(SurvivorSelectionOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(SurvivorSelectionOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_sselect_is_abstract():
    generic_assert_abstract_method(SurvivorSelectionOperator, "sselect")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_SurvivorSelectionOperator_is_concrete():
    generic_assert_concrete_function(is_SurvivorSelectionOperator)

def test_check_is_SurvivorSelectionOperator_is_concrete():
    generic_assert_concrete_function(check_is_SurvivorSelectionOperator)

def test_cond_check_is_SurvivorSelectionOperator_is_concrete():
    generic_assert_concrete_function(cond_check_is_SurvivorSelectionOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_SurvivorSelectionOperator(operator):
    assert is_SurvivorSelectionOperator(operator)

def test_check_is_SurvivorSelectionOperator(operator):
    with not_raises(TypeError):
        check_is_SurvivorSelectionOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_SurvivorSelectionOperator(None, "operator")