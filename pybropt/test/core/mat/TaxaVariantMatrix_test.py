import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import TaxaVariantMatrix
from pybropt.core.mat import is_TaxaVariantMatrix
from pybropt.core.mat import check_is_TaxaVariantMatrix
from pybropt.core.mat import cond_check_is_TaxaVariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield TaxaVariantMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(TaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(TaxaVariantMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_TaxaVariantMatrix_is_concrete():
    generic_assert_concrete_function(is_TaxaVariantMatrix)

def test_is_TaxaVariantMatrix(mat):
    assert is_TaxaVariantMatrix(mat)

def test_check_is_TaxaVariantMatrix_is_concrete():
    generic_assert_concrete_function(check_is_TaxaVariantMatrix)

def test_check_is_TaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaVariantMatrix(None, "mat")

def test_cond_check_is_TaxaVariantMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_TaxaVariantMatrix)
