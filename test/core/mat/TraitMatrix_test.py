import inspect
import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.TraitMatrix import TraitMatrix
from pybrops.core.mat.TraitMatrix import check_is_TraitMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyTraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(TraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_trait_is_abstract():
    assert_property_isabstract(TraitMatrix, "trait")

def test_ntrait_is_abstract():
    assert_property_isabstract(TraitMatrix, "ntrait")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "adjoin_trait")

def test_delete_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "delete_trait")

def test_insert_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "insert_trait")

def test_select_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "select_trait")

def test_concat_trait_is_abstract():
    assert_classmethod_isabstract(TraitMatrix, "concat_trait")

def test_append_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "append_trait")

def test_remove_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "remove_trait")

def test_incorp_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "incorp_trait")

def test_lexsort_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "lexsort_trait")

def test_reorder_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "reorder_trait")

def test_sort_trait_is_abstract():
    assert_method_isabstract(TraitMatrix, "sort_trait")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_TraitMatrix_is_concrete():
    assert_function_isconcrete(check_is_TraitMatrix)

def test_check_is_TraitMatrix(mat):
    with not_raises(TypeError):
        check_is_TraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TraitMatrix(None, "mat")
