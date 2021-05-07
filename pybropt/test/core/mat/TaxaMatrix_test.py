import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import TaxaMatrix
from pybropt.core.mat import is_TaxaMatrix
from pybropt.core.mat import check_is_TaxaMatrix
from pybropt.core.mat import cond_check_is_TaxaMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield TaxaMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(TaxaMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(TaxaMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_taxa_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "taxa")

def test_ntaxa_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "ntaxa")

def test_taxa_grp_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "taxa_grp")

def test_taxa_grp_name_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "taxa_grp_name")

def test_taxa_grp_stix_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "taxa_grp_stix")

def test_taxa_grp_spix_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "taxa_grp_spix")

def test_taxa_grp_len_is_abstract():
    generic_assert_abstract_property(TaxaMatrix, "taxa_grp_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract(mat):
    generic_assert_abstract_method(mat, "lexsort")

def test_adjoin_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "adjoin_taxa")

def test_delete_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "delete_taxa")

def test_insert_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "insert_taxa")

def test_select_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "select_taxa")

def test_concat_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "concat_taxa")

def test_append_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "append_taxa")

def test_remove_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "remove_taxa")

def test_incorp_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "incorp_taxa")

def test_lexsort_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "lexsort_taxa")

def test_reorder_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "reorder_taxa")

def test_sort_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "sort_taxa")

def test_group_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "group_taxa")

def test_is_grouped_taxa_is_abstract(mat):
    generic_assert_abstract_method(mat, "is_grouped_taxa")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_TaxaMatrix_is_concrete():
    generic_assert_concrete_function(is_TaxaMatrix)

def test_is_TaxaMatrix(mat):
    assert is_TaxaMatrix(mat)

def test_check_is_TaxaMatrix_is_concrete():
    generic_assert_concrete_function(check_is_TaxaMatrix)

def test_check_is_TaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaMatrix(None, "mat")

def test_cond_check_is_TaxaMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_TaxaMatrix)
