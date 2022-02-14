import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.core.mat.MutableMatrix import MutableMatrix
from pybrops.core.mat.MutableMatrix import is_MutableMatrix
from pybrops.core.mat.MutableMatrix import check_is_MutableMatrix
from pybrops.core.mat.MutableMatrix import cond_check_is_MutableMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield MutableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(MutableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(MutableMatrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract(mat):
    generic_assert_abstract_method(mat, "append")

def test_remove_is_abstract(mat):
    generic_assert_abstract_method(mat, "remove")

def test_incorp_is_abstract(mat):
    generic_assert_abstract_method(mat, "incorp")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_MutableMatrix_is_concrete():
    generic_assert_concrete_function(is_MutableMatrix)

def test_is_MutableMatrix(mat):
    assert is_MutableMatrix(mat)

def test_check_is_MutableMatrix_is_concrete():
    generic_assert_concrete_function(check_is_MutableMatrix)

def test_check_is_MutableMatrix(mat):
    with not_raises(TypeError):
        check_is_MutableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_MutableMatrix(None, "mat")

def test_cond_check_is_MutableMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_MutableMatrix)
