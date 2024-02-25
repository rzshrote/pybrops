import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.arch.ImmigrationOperator import ImmigrationOperator
from pybrops.breed.arch.ImmigrationOperator import check_is_ImmigrationOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyImmigrationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(ImmigrationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(ImmigrationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(ImmigrationOperator, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_ImmigrationOperator_is_concrete():
    assert_function_isconcrete(check_is_ImmigrationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_ImmigrationOperator(arch):
    with not_raises(TypeError):
        check_is_ImmigrationOperator(arch, "arch")
    with pytest.raises(TypeError):
        check_is_ImmigrationOperator(None, "arch")
