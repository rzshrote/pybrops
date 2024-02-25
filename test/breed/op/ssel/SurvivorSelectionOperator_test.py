import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import check_is_SurvivorSelectionOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield DummySurvivorSelectionOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(SurvivorSelectionOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(SurvivorSelectionOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_sselect_is_abstract():
    assert_method_isabstract(SurvivorSelectionOperator, "sselect")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_SurvivorSelectionOperator_is_concrete():
    assert_function_isconcrete(check_is_SurvivorSelectionOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SurvivorSelectionOperator(operator):
    with not_raises(TypeError):
        check_is_SurvivorSelectionOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_SurvivorSelectionOperator(None, "operator")
