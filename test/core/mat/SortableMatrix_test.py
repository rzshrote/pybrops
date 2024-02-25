import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.SortableMatrix import SortableMatrix
from pybrops.core.mat.SortableMatrix import check_is_SortableMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummySortableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(SortableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract():
    assert_method_isabstract(SortableMatrix, "lexsort")

def test_reorder_is_abstract():
    assert_method_isabstract(SortableMatrix, "reorder")

def test_sort_is_abstract():
    assert_method_isabstract(SortableMatrix, "sort")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SortableMatrix_is_concrete():
    assert_function_isconcrete(check_is_SortableMatrix)

def test_check_is_SortableMatrix(mat):
    with not_raises(TypeError):
        check_is_SortableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SortableMatrix(None, "mat")
