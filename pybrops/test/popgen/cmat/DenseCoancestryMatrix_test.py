import pytest
import numpy

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.cmat.DenseCoancestryMatrix import is_DenseCoancestryMatrix
from pybrops.popgen.cmat.DenseCoancestryMatrix import check_is_DenseCoancestryMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def X_mat_int8():
    X = numpy.array(
       [[-1,  1,  0,  0, -1, -1,  0, -1,  0, -1, -1,  1, -1,  0,  1, -1],
        [-1,  1, -1, -1,  0,  1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1],
        [ 1,  1,  1,  0,  0,  1,  1, -1, -1,  0,  0,  1, -1, -1, -1, -1],
        [-1,  1,  1, -1,  1,  1,  1,  0,  1,  1,  0,  1,  1, -1, -1,  1],
        [ 1,  1,  1, -1,  0,  0,  0, -1, -1,  1, -1,  0,  1,  0, -1, -1],
        [-1, -1,  1,  0, -1,  0,  1, -1, -1,  0, -1,  1,  0, -1, -1,  1],
        [ 0,  0,  0,  1, -1,  0,  1, -1,  1,  0, -1,  1, -1,  0, -1,  0],
        [-1,  1,  0, -1,  0,  1,  1,  0,  0, -1, -1,  1, -1, -1,  0, -1]], 
        dtype = "int8"
    )
    return X

@pytest.fixture
def A_mat_float64(X_mat_int8):
    A = ((1.0/X_mat_int8.shape[1]) * (X_mat_int8.dot(X_mat_int8.T))) + 1.0
    yield A

@pytest.fixture
def cmat(A_mat_float64):
    yield DenseCoancestryMatrix(A_mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseCoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "__init__")

def test_coancestry_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "coancestry")

def test_kinship_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "kinship")

def test_is_positive_semidefinite_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "is_positive_semidefinite")

def test_apply_jitter_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "apply_jitter")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_gmat_is_abstract():
    assert_abstract_method(DenseCoancestryMatrix, "from_gmat")

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################
def test_mat_fget(cmat, A_mat_float64):
    assert numpy.all(cmat == A_mat_float64)

def test_mat_fset_TypeError(cmat, A_mat_float64):
    with pytest.raises(TypeError):
        cmat.mat = list(A_mat_float64.flatten())

def test_mat_fset_ValueError(cmat, A_mat_float64):
    with pytest.raises(ValueError):
        cmat.mat = A_mat_float64.flatten()

def test_mat_fset(cmat, A_mat_float64):
    cmat.mat = A_mat_float64
    assert numpy.all(cmat.mat == A_mat_float64)

def test_mat_fdel(cmat, A_mat_float64):
    del cmat.mat
    with pytest.raises(AttributeError):
        cmat.mat

################################################################################
###################### Test concrete method functionality ######################
################################################################################

# matrix conversion
def test_mat_asformat_TypeError(cmat):
    with pytest.raises(TypeError):
        K = cmat.mat_asformat([])

def test_mat_asformat_ValueError(cmat):
    with pytest.raises(ValueError):
        K = cmat.mat_asformat("unknown")

def test_mat_asformat(cmat, A_mat_float64):
    A = cmat.mat_asformat("coancestry")
    assert numpy.all(A == A_mat_float64)
    K = cmat.mat_asformat("kinship")
    assert numpy.all(K == (0.5 * A_mat_float64))

# coancestry function
def test_coancestry_2tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.coancestry(i,j) == A_mat_float64[i,j]

def test_coancestry_1tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i) == A_mat_float64[i])

def test_coancestry_row_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i,slice(None)) == A_mat_float64[i,:])

def test_coancestry_col_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(slice(None),i) == A_mat_float64[:,i])

def test_coancestry_list_tuple(cmat, A_mat_float64):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.coancestry(a,b) == A_mat_float64[a,b])

# kinship function
def test_kinship_2tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.kinship(i,j) == (0.5 * A_mat_float64[i,j])

def test_kinship_1tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(i) == (0.5 * A_mat_float64[i]))

def test_kinship_row_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(i,slice(None)) == (0.5 * A_mat_float64[i,:]))

def test_kinship_col_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(slice(None),i) == (0.5 * A_mat_float64[:,i]))

def test_kinship_list_tuple(cmat, A_mat_float64):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.kinship(a,b) == (0.5 * A_mat_float64[a,b]))

# postitive definite checks
def test_is_positive_semidefinite(cmat):
    assert cmat.is_positive_semidefinite()

def test_is_positive_semidefinite_eigvaltol(cmat):
    assert cmat.is_positive_semidefinite(-1.0)

# jitter function checks
def test_apply_jitter(cmat):
    assert cmat.apply_jitter()

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_DenseCoancestryMatrix_is_concrete():
    assert_concrete_function(is_DenseCoancestryMatrix)

def test_check_is_DenseCoancestryMatrix_is_concrete():
    assert_concrete_function(check_is_DenseCoancestryMatrix)

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
