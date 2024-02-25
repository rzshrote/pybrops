import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.MutableMatrix import MutableMatrix
from pybrops.core.mat.MutableMatrix import check_is_MutableMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyMutableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(MutableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract():
    assert_method_isabstract(MutableMatrix, "append")

def test_remove_is_abstract():
    assert_method_isabstract(MutableMatrix, "remove")

def test_incorp_is_abstract():
    assert_method_isabstract(MutableMatrix, "incorp")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_MutableMatrix_is_concrete():
    assert_function_isconcrete(check_is_MutableMatrix)

def test_check_is_MutableMatrix(mat):
    with not_raises(TypeError):
        check_is_MutableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_MutableMatrix(None, "mat")
