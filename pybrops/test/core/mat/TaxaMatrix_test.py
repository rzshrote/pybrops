import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.TaxaMatrix import check_is_TaxaMatrix
from pybrops.test.core.mat.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyTaxaMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(TaxaMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(TaxaMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_taxa_is_abstract():
    assert_abstract_property(TaxaMatrix, "taxa")

def test_ntaxa_is_abstract():
    assert_abstract_property(TaxaMatrix, "ntaxa")

def test_taxa_grp_is_abstract():
    assert_abstract_property(TaxaMatrix, "taxa_grp")

def test_taxa_grp_name_is_abstract():
    assert_abstract_property(TaxaMatrix, "taxa_grp_name")

def test_taxa_grp_stix_is_abstract():
    assert_abstract_property(TaxaMatrix, "taxa_grp_stix")

def test_taxa_grp_spix_is_abstract():
    assert_abstract_property(TaxaMatrix, "taxa_grp_spix")

def test_taxa_grp_len_is_abstract():
    assert_abstract_property(TaxaMatrix, "taxa_grp_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_lexsort_is_abstract():
    assert_abstract_method(TaxaMatrix, "lexsort")

def test_adjoin_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "adjoin_taxa")

def test_delete_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "delete_taxa")

def test_insert_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "insert_taxa")

def test_select_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "select_taxa")

def test_concat_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "concat_taxa")

def test_append_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "append_taxa")

def test_remove_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "remove_taxa")

def test_incorp_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "incorp_taxa")

def test_lexsort_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "lexsort_taxa")

def test_reorder_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "reorder_taxa")

def test_sort_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "sort_taxa")

def test_group_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "group_taxa")

def test_is_grouped_taxa_is_abstract():
    assert_abstract_method(TaxaMatrix, "is_grouped_taxa")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_TaxaMatrix_is_concrete():
    assert_concrete_function(check_is_TaxaMatrix)

def test_check_is_TaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaMatrix(None, "mat")
