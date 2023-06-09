import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.mate.MatingOperator import check_is_MatingOperator

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def operator():
    yield MatingOperator()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(MatingOperator)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(MatingOperator, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_mate_is_abstract():
    assert_abstract_method(MatingOperator, "mate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_MatingOperator_is_concrete():
    assert_concrete_function(check_is_MatingOperator)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_MatingOperator(operator):
    with not_raises(TypeError):
        check_is_MatingOperator(operator, "operator")
    with pytest.raises(TypeError):
        check_is_MatingOperator(None, "operator")
