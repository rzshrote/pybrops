import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.core.mat.SortableMatrix import SortableMatrix
from pybrops.core.mat.SortableMatrix import is_SortableMatrix
from pybrops.core.mat.SortableMatrix import check_is_SortableMatrix

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
    assert_docstring(SortableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SortableMatrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract(mat):
    assert_abstract_method(mat, "lexsort")

def test_reorder_is_abstract(mat):
    assert_abstract_method(mat, "reorder")

def test_sort_is_abstract(mat):
    assert_abstract_method(mat, "sort")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_SortableMatrix_is_concrete():
    assert_concrete_function(is_SortableMatrix)

def test_is_SortableMatrix(mat):
    assert is_SortableMatrix(mat)

def test_check_is_SortableMatrix_is_concrete():
    assert_concrete_function(check_is_SortableMatrix)

def test_check_is_SortableMatrix(mat):
    with not_raises(TypeError):
        check_is_SortableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SortableMatrix(None, "mat")
