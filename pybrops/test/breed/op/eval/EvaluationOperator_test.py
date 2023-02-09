import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator
from pybrops.breed.op.eval.EvaluationOperator import is_EvaluationOperator
from pybrops.breed.op.eval.EvaluationOperator import check_is_EvaluationOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield EvaluationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(EvaluationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(EvaluationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_evaluate_is_abstract():
    generic_assert_abstract_method(EvaluationOperator, "evaluate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_EvaluationOperator_is_concrete():
    generic_assert_concrete_function(is_EvaluationOperator)

def test_check_is_EvaluationOperator_is_concrete():
    generic_assert_concrete_function(check_is_EvaluationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_EvaluationOperator(operator):
    assert is_EvaluationOperator(operator)

def test_check_is_EvaluationOperator(operator):
    with not_raises(TypeError):
        check_is_EvaluationOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_EvaluationOperator(None, "operator")
