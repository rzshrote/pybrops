import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybropt.core.mat.PhasedTaxaVariantMatrix import is_PhasedTaxaVariantMatrix
from pybropt.core.mat.PhasedTaxaVariantMatrix import check_is_PhasedTaxaVariantMatrix
from pybropt.core.mat.PhasedTaxaVariantMatrix import cond_check_is_PhasedTaxaVariantMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield PhasedTaxaVariantMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(PhasedTaxaVariantMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(PhasedTaxaVariantMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhasedTaxaVariantMatrix_is_concrete():
    generic_assert_concrete_function(is_PhasedTaxaVariantMatrix)

def test_is_PhasedTaxaVariantMatrix(mat):
    assert is_PhasedTaxaVariantMatrix(mat)

def test_check_is_PhasedTaxaVariantMatrix_is_concrete():
    generic_assert_concrete_function(check_is_PhasedTaxaVariantMatrix)

def test_check_is_PhasedTaxaVariantMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedTaxaVariantMatrix(None, "mat")

def test_cond_check_is_PhasedTaxaVariantMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_PhasedTaxaVariantMatrix)
