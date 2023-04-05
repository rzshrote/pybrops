import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.core.mat.VariantMatrix import VariantMatrix
from pybrops.core.mat.VariantMatrix import is_VariantMatrix
from pybrops.core.mat.VariantMatrix import check_is_VariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield VariantMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(VariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(VariantMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_vrnt_chrgrp_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_chrgrp")

def test_vrnt_phypos_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_phypos")

def test_vrnt_name_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_name")

def test_vrnt_genpos_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_genpos")

def test_vrnt_xoprob_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_xoprob")

def test_vrnt_hapgrp_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_hapgrp")

def test_vrnt_mask_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_mask")

def test_vrnt_chrgrp_name_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_chrgrp_name")

def test_vrnt_chrgrp_stix_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_chrgrp_stix")

def test_vrnt_chrgrp_spix_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_chrgrp_spix")

def test_vrnt_chrgrp_len_is_abstract():
    assert_abstract_property(VariantMatrix, "vrnt_chrgrp_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "adjoin_vrnt")

def test_delete_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "delete_vrnt")

def test_insert_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "insert_vrnt")

def test_select_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "select_vrnt")

def test_concat_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "concat_vrnt")

def test_append_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "append_vrnt")

def test_remove_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "remove_vrnt")

def test_incorp_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "incorp_vrnt")

def test_lexsort_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "lexsort_vrnt")

def test_reorder_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "reorder_vrnt")

def test_sort_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "sort_vrnt")

def test_group_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "group_vrnt")

def test_is_grouped_vrnt_is_abstract(mat):
    assert_abstract_method(mat, "is_grouped_vrnt")

def test_interp_genpos_is_abstract(mat):
    assert_abstract_method(mat, "interp_genpos")

def test_interp_xoprob_is_abstract(mat):
    assert_abstract_method(mat, "interp_xoprob")

def test_assign_hapgrp_is_abstract(mat):
    assert_abstract_method(mat, "assign_hapgrp")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_VariantMatrix_is_concrete():
    assert_concrete_function(is_VariantMatrix)

def test_is_VariantMatrix(mat):
    assert is_VariantMatrix(mat)

def test_check_is_VariantMatrix_is_concrete():
    assert_concrete_function(check_is_VariantMatrix)

def test_check_is_VariantMatrix(mat):
    with not_raises(TypeError):
        check_is_VariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_VariantMatrix(None, "mat")
