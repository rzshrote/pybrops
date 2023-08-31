import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.MutableMatrix import MutableMatrix
from pybrops.core.mat.MutableMatrix import check_is_MutableMatrix
from pybrops.test.core.mat.common_fixtures import *

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
    assert_docstring(MutableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(MutableMatrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract():
    assert_abstract_method(MutableMatrix, "append")

def test_remove_is_abstract():
    assert_abstract_method(MutableMatrix, "remove")

def test_incorp_is_abstract():
    assert_abstract_method(MutableMatrix, "incorp")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_MutableMatrix_is_concrete():
    assert_concrete_function(check_is_MutableMatrix)

def test_check_is_MutableMatrix(mat):
    with not_raises(TypeError):
        check_is_MutableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_MutableMatrix(None, "mat")
