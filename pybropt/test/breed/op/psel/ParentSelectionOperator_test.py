import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.op.psel import ParentSelectionOperator
from pybropt.breed.op.psel import is_ParentSelectionOperator
from pybropt.breed.op.psel import check_is_ParentSelectionOperator
from pybropt.breed.op.psel import cond_check_is_ParentSelectionOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield ParentSelectionOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(ParentSelectionOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(ParentSelectionOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_pselect_is_abstract():
    generic_assert_abstract_method(ParentSelectionOperator, "pselect")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_ParentSelectionOperator_is_concrete():
    generic_assert_concrete_function(is_ParentSelectionOperator)

def test_check_is_ParentSelectionOperator_is_concrete():
    generic_assert_concrete_function(check_is_ParentSelectionOperator)

def test_cond_check_is_ParentSelectionOperator_is_concrete():
    generic_assert_concrete_function(cond_check_is_ParentSelectionOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_ParentSelectionOperator(operator):
    assert is_ParentSelectionOperator(operator)

def test_check_is_ParentSelectionOperator(operator):
    with not_raises(TypeError):
        check_is_ParentSelectionOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_ParentSelectionOperator(None, "operator")