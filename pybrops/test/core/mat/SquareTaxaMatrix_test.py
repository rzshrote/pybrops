import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.SquareTaxaMatrix import is_SquareTaxaMatrix
from pybrops.core.mat.SquareTaxaMatrix import check_is_SquareTaxaMatrix
from pybrops.core.mat.SquareTaxaMatrix import cond_check_is_SquareTaxaMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield SquareTaxaMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(SquareTaxaMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(SquareTaxaMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_SquareTaxaMatrix_is_concrete():
    generic_assert_concrete_function(is_SquareTaxaMatrix)

def test_is_SquareTaxaMatrix(mat):
    assert is_SquareTaxaMatrix(mat)

def test_check_is_SquareTaxaMatrix_is_concrete():
    generic_assert_concrete_function(check_is_SquareTaxaMatrix)

def test_check_is_SquareTaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareTaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareTaxaMatrix(None, "mat")

def test_cond_check_is_SquareTaxaMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_SquareTaxaMatrix)
