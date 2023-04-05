import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.core.mat.TraitMatrix import TraitMatrix
from pybrops.core.mat.TraitMatrix import is_TraitMatrix
from pybrops.core.mat.TraitMatrix import check_is_TraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield TraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(TraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(TraitMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_trait_is_abstract():
    assert_abstract_property(TraitMatrix, "trait")

def test_ntrait_is_abstract():
    assert_abstract_property(TraitMatrix, "ntrait")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_trait_is_abstract(mat):
    assert_abstract_method(mat, "adjoin_trait")

def test_delete_trait_is_abstract(mat):
    assert_abstract_method(mat, "delete_trait")

def test_insert_trait_is_abstract(mat):
    assert_abstract_method(mat, "insert_trait")

def test_select_trait_is_abstract(mat):
    assert_abstract_method(mat, "select_trait")

def test_concat_trait_is_abstract(mat):
    assert_abstract_method(mat, "concat_trait")

def test_append_trait_is_abstract(mat):
    assert_abstract_method(mat, "append_trait")

def test_remove_trait_is_abstract(mat):
    assert_abstract_method(mat, "remove_trait")

def test_incorp_trait_is_abstract(mat):
    assert_abstract_method(mat, "incorp_trait")

def test_lexsort_trait_is_abstract(mat):
    assert_abstract_method(mat, "lexsort_trait")

def test_reorder_trait_is_abstract(mat):
    assert_abstract_method(mat, "reorder_trait")

def test_sort_trait_is_abstract(mat):
    assert_abstract_method(mat, "sort_trait")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_TraitMatrix_is_concrete():
    assert_concrete_function(is_TraitMatrix)

def test_is_TraitMatrix(mat):
    assert is_TraitMatrix(mat)

def test_check_is_TraitMatrix_is_concrete():
    assert_concrete_function(check_is_TraitMatrix)

def test_check_is_TraitMatrix(mat):
    with not_raises(TypeError):
        check_is_TraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TraitMatrix(None, "mat")
