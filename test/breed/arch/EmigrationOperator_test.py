import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.arch.EmigrationOperator import EmigrationOperator
from pybrops.breed.arch.EmigrationOperator import check_is_EmigrationOperator
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyEmigrationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(EmigrationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(EmigrationOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(EmigrationOperator, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_EmigrationOperator_is_concrete():
    assert_function_isconcrete(check_is_EmigrationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_EmigrationOperator(arch):
    with not_raises(TypeError):
        check_is_EmigrationOperator(arch, "arch")
    with pytest.raises(TypeError):
        check_is_EmigrationOperator(None, "arch")
