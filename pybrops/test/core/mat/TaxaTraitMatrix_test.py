import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.core.mat.TaxaTraitMatrix import TaxaTraitMatrix
from pybrops.core.mat.TaxaTraitMatrix import is_TaxaTraitMatrix
from pybrops.core.mat.TaxaTraitMatrix import check_is_TaxaTraitMatrix
from pybrops.core.mat.TaxaTraitMatrix import cond_check_is_TaxaTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield TaxaTraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(TaxaTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(TaxaTraitMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_TaxaTraitMatrix_is_concrete():
    generic_assert_concrete_function(is_TaxaTraitMatrix)

def test_is_TaxaTraitMatrix(mat):
    assert is_TaxaTraitMatrix(mat)

def test_check_is_TaxaTraitMatrix_is_concrete():
    generic_assert_concrete_function(check_is_TaxaTraitMatrix)

def test_check_is_TaxaTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_TaxaTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_TaxaTraitMatrix(None, "mat")

def test_cond_check_is_TaxaTraitMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_TaxaTraitMatrix)
