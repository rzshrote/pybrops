import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.arch.EmigrationOperator import EmigrationOperator
from pybrops.breed.arch.EmigrationOperator import check_is_EmigrationOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield EmigrationOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(EmigrationOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(EmigrationOperator, "__init__")

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
    assert_concrete_function(check_is_EmigrationOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_EmigrationOperator(arch):
    with not_raises(TypeError):
        check_is_EmigrationOperator(arch, "arch")
    with pytest.raises(TypeError):
        check_is_EmigrationOperator(None, "arch")
