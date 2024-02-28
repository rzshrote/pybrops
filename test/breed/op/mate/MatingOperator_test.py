import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.mate.MatingOperator import check_is_MatingOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield DummyMatingOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(MatingOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mate_is_abstract():
    assert_method_isabstract(MatingOperator, "mate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_MatingOperator_is_concrete():
    assert_function_isconcrete(check_is_MatingOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_MatingOperator(operator):
    with not_raises(TypeError):
        check_is_MatingOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_MatingOperator(None, "operator")
