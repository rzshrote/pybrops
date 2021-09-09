import pytest
import numpy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import DenseMatrix
from pybropt.core.mat import is_DenseMatrix
from pybropt.core.mat import check_is_DenseMatrix
from pybropt.core.mat import cond_check_is_DenseMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [3.3, 9.2, 5.6],
        [8.7, 3.7, 4.1],
        [9.0, 4.7, 3.8]
    ])
    yield a

@pytest.fixture
def mat(mat_float64):
    yield DenseMatrix(mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DenseMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseMatrix, "__init__")

def test_add_is_concrete():
    generic_assert_concrete_method(DenseMatrix, "__add__")

def test_add_is_concrete():
    generic_assert_concrete_method(DenseMatrix, "__sub__")

def test_add_is_concrete():
    generic_assert_concrete_method(DenseMatrix, "__mul__")

def test_add_is_concrete():
    generic_assert_concrete_method(DenseMatrix, "__matmul__")

################################################################################
############################ Test Class Properties #############################
################################################################################
def test_mat_fget(mat, mat_float64):
    assert numpy.all(mat.mat == mat_float64)

################################################################################
######################### Test special class functions #########################
################################################################################
def test_add(mat, mat_float64):
    assert numpy.all((mat + mat_float64) == (mat_float64 + mat_float64))
    assert numpy.all((mat + mat) == (mat_float64 + mat_float64))

def test_sub(mat, mat_float64):
    assert numpy.all((mat - mat_float64) == (mat_float64 - mat_float64))
    assert numpy.all((mat - mat) == (mat_float64 - mat_float64))

def test_mul(mat, mat_float64):
    assert numpy.all((mat * mat_float64) == (mat_float64 * mat_float64))
    assert numpy.all((mat * mat) == (mat_float64 * mat_float64))

def test_matmul(mat, mat_float64):
    assert numpy.all((mat @ mat_float64) == (mat_float64 @ mat_float64))
    assert numpy.all((mat @ mat) == (mat_float64 @ mat_float64))

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseMatrix)

def test_is_DenseMatrix(mat):
    assert is_DenseMatrix(mat)

def test_check_is_DenseMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseMatrix)

def test_check_is_DenseMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseMatrix(None, "mat")

def test_cond_check_is_DenseMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseMatrix)
