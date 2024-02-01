import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator
from pybrops.breed.op.eval.EvaluationOperator import check_is_EvaluationOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield DummyEvaluationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(EvaluationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(EvaluationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_evaluate_is_abstract():
    assert_abstract_method(EvaluationOperator, "evaluate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_EvaluationOperator_is_concrete():
    assert_concrete_function(check_is_EvaluationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_EvaluationOperator(operator):
    with not_raises(TypeError):
        check_is_EvaluationOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_EvaluationOperator(None, "operator")
