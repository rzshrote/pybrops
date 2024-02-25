import inspect
import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.VariantMatrix import VariantMatrix
from pybrops.core.mat.VariantMatrix import check_is_VariantMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyVariantMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(VariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_vrnt_chrgrp_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_chrgrp")

def test_vrnt_phypos_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_phypos")

def test_vrnt_name_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_name")

def test_vrnt_genpos_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_genpos")

def test_vrnt_xoprob_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_xoprob")

def test_vrnt_hapgrp_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_hapgrp")

def test_vrnt_mask_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_mask")

def test_vrnt_chrgrp_name_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_chrgrp_name")

def test_vrnt_chrgrp_stix_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_chrgrp_stix")

def test_vrnt_chrgrp_spix_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_chrgrp_spix")

def test_vrnt_chrgrp_len_is_abstract():
    assert_property_isabstract(VariantMatrix, "vrnt_chrgrp_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "adjoin_vrnt")

def test_delete_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "delete_vrnt")

def test_insert_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "insert_vrnt")

def test_select_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "select_vrnt")

def test_concat_vrnt_is_abstract():
    assert_classmethod_isabstract(VariantMatrix, "concat_vrnt")

def test_append_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "append_vrnt")

def test_remove_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "remove_vrnt")

def test_incorp_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "incorp_vrnt")

def test_lexsort_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "lexsort_vrnt")

def test_reorder_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "reorder_vrnt")

def test_sort_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "sort_vrnt")

def test_group_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "group_vrnt")

def test_is_grouped_vrnt_is_abstract():
    assert_method_isabstract(VariantMatrix, "is_grouped_vrnt")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_VariantMatrix_is_concrete():
    assert_function_isconcrete(check_is_VariantMatrix)

def test_check_is_VariantMatrix(mat):
    with not_raises(TypeError):
        check_is_VariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_VariantMatrix(None, "mat")
