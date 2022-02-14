import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.SquareMatrix import is_SquareMatrix
from pybrops.core.mat.SquareMatrix import check_is_SquareMatrix
from pybrops.core.mat.SquareMatrix import cond_check_is_SquareMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield SquareMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(SquareMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(SquareMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_mat_is_abstract():
    generic_assert_abstract_property(SquareMatrix, "nsquare")

def test_mat_is_abstract():
    generic_assert_abstract_property(SquareMatrix, "square_axes")

def test_mat_is_abstract():
    generic_assert_abstract_property(SquareMatrix, "square_axes_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract(mat):
    generic_assert_abstract_method(mat, "is_square")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_SquareMatrix_is_concrete():
    generic_assert_concrete_function(is_SquareMatrix)

def test_is_SquareMatrix(mat):
    assert is_SquareMatrix(mat)

def test_check_is_SquareMatrix_is_concrete():
    generic_assert_concrete_function(check_is_SquareMatrix)

def test_check_is_SquareMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareMatrix(None, "mat")

def test_cond_check_is_SquareMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_SquareMatrix)
