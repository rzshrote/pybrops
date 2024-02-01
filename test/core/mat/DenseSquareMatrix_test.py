import pytest
import numpy
import copy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.DenseSquareMatrix import DenseSquareMatrix
from pybrops.core.mat.DenseSquareMatrix import check_is_DenseSquareMatrix

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
    out = DenseSquareMatrix(mat_float64)
    yield out

@pytest.fixture
def mat_rectangle():
    a = numpy.float64([
        [3.3, 5.6],
        [8.7, 4.1],
        [9.0, 3.8]
    ])
    out = DenseSquareMatrix(a)
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseSquareMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################## Matrix Data Properties ##################
def test_mat_fget(mat, mat_float64):
    assert numpy.all(mat.mat == mat_float64)

def test_mat_fset(mat, mat_float64):
    mat.mat = mat_float64.T
    assert numpy.all(mat.mat == mat_float64.T)

def test_mat_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.mat

################ Square Metadata Properites ################
def test_nsquare_fget(mat, mat_float64):
    assert mat.nsquare == len(mat.square_axes)

def test_nsquare_fset(mat, mat_float64):
    with pytest.raises(AttributeError):
        mat.nsquare = len(mat_float64)

def test_nsquare_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.nsquare

def test_square_axes_fget(mat):
    assert mat.square_axes == (0,1)

def test_square_axes_fset(mat):
    with pytest.raises(AttributeError):
        mat.square_axes = (3,4)

def test_square_axes_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.square_axes

def test_square_axes_len_fget(mat, mat_float64):
    assert numpy.all(mat.square_axes_len == mat_float64.shape)

def test_square_axes_len_fset(mat, mat_float64):
    with pytest.raises(AttributeError):
        mat.square_axes_len = (5,5)

def test_square_axes_len_fdel(mat):
    with pytest.raises(AttributeError):
        del mat.square_axes_len

##################### Fill data lookup #####################
def test_fill_value(mat):
    assert isinstance(mat._fill_value, dict)

################################################################################
###################### Test concrete method functionality ######################
################################################################################

###################### Square Methods ######################
def test_is_square(mat, mat_rectangle):
    assert mat.is_square()
    assert not mat_rectangle.is_square()

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseSquareMatrix_is_concrete():
    assert_concrete_function(check_is_DenseSquareMatrix)

def test_check_is_DenseSquareMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseSquareMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseSquareMatrix(None, "mat")
