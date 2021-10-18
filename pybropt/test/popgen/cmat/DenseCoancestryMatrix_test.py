import pytest
import numpy

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.popgen.cmat import DenseCoancestryMatrix
from pybropt.popgen.cmat import is_DenseCoancestryMatrix
from pybropt.popgen.cmat import check_is_DenseCoancestryMatrix
from pybropt.popgen.cmat import cond_check_is_DenseCoancestryMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [0.57, 0.18, 0.54, 0.45, 0.09, 0.6 , 0.63, 0.4 ],
        [0.18, 0.77, 0.8 , 0.53, 0.89, 0.34, 0.41, 0.71],
        [0.54, 0.8 , 0.46, 0.52, 0.65, 0.85, 0.21, 0.68],
        [0.45, 0.53, 0.52, 0.73, 0.54, 0.64, 0.41, 0.64],
        [0.09, 0.89, 0.65, 0.54, 0.3 , 0.62, 0.51, 0.74],
        [0.6 , 0.34, 0.85, 0.64, 0.62, 0.96, 0.53, 0.45],
        [0.63, 0.41, 0.21, 0.41, 0.51, 0.53, 0.45, 0.38],
        [0.4 , 0.71, 0.68, 0.64, 0.74, 0.45, 0.38, 0.65]
    ])
    yield a

@pytest.fixture
def cmat(mat_float64):
    yield DenseCoancestryMatrix(mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(DenseCoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseCoancestryMatrix, "__init__")

def test_coancestry_is_concrete():
    generic_assert_concrete_method(DenseCoancestryMatrix, "coancestry")

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################
def test_mat_fget(cmat, mat_float64):
    assert numpy.all(cmat == mat_float64)

def test_mat_fset_TypeError(cmat, mat_float64):
    with pytest.raises(TypeError):
        cmat.mat = list(mat_float64.flatten())

def test_mat_fset_ValueError(cmat, mat_float64):
    with pytest.raises(ValueError):
        cmat.mat = mat_float64.flatten()

def test_mat_fset(cmat, mat_float64):
    cmat.mat = mat_float64
    assert numpy.all(cmat.mat == mat_float64)

def test_mat_fdel(cmat, mat_float64):
    del cmat.mat
    with pytest.raises(AttributeError):
        cmat.mat

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_coancestry_2tuple(cmat, mat_float64):
    n = mat_float64.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.coancestry(i,j) == mat_float64[i,j]

def test_coancestry_1tuple(cmat, mat_float64):
    n = mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i) == mat_float64[i])

def test_coancestry_row_slice(cmat, mat_float64):
    n = mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i,slice(None)) == mat_float64[i,:])

def test_coancestry_col_slice(cmat, mat_float64):
    n = mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(slice(None),i) == mat_float64[:,i])

def test_coancestry_list_tuple(cmat, mat_float64):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.coancestry(a,b) == mat_float64[a,b])

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DenseCoancestryMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseCoancestryMatrix)

def test_check_is_DenseCoancestryMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseCoancestryMatrix)

def test_cond_check_is_DenseCoancestryMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseCoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseCoancestryMatrix(cmat):
    assert is_DenseCoancestryMatrix(cmat)

def test_check_is_DenseCoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_DenseCoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_DenseCoancestryMatrix(None, "cmat")
