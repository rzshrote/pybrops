import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

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
    assert_class_documentation(EvaluationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(EvaluationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_evaluate_is_abstract():
    assert_method_isabstract(EvaluationOperator, "evaluate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_EvaluationOperator_is_concrete():
    assert_function_isconcrete(check_is_EvaluationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_EvaluationOperator(operator):
    with not_raises(TypeError):
        check_is_EvaluationOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_EvaluationOperator(None, "operator")
