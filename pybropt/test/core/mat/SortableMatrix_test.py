import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat.SortableMatrix import SortableMatrix
from pybropt.core.mat.SortableMatrix import is_SortableMatrix
from pybropt.core.mat.SortableMatrix import check_is_SortableMatrix
from pybropt.core.mat.SortableMatrix import cond_check_is_SortableMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield SortableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(SortableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(SortableMatrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract(mat):
    generic_assert_abstract_method(mat, "lexsort")

def test_reorder_is_abstract(mat):
    generic_assert_abstract_method(mat, "reorder")

def test_sort_is_abstract(mat):
    generic_assert_abstract_method(mat, "sort")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_SortableMatrix_is_concrete():
    generic_assert_concrete_function(is_SortableMatrix)

def test_is_SortableMatrix(mat):
    assert is_SortableMatrix(mat)

def test_check_is_SortableMatrix_is_concrete():
    generic_assert_concrete_function(check_is_SortableMatrix)

def test_check_is_SortableMatrix(mat):
    with not_raises(TypeError):
        check_is_SortableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SortableMatrix(None, "mat")

def test_cond_check_is_SortableMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_SortableMatrix)
