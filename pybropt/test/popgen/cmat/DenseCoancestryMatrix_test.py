import numpy
import pytest

from pybropt.popgen.cmat import DenseCoancestryMatrix
from pybropt.popgen.cmat import is_DenseCoancestryMatrix

from pybropt.test import generic_test_ndarray_add
from pybropt.test import generic_test_ndarray_sub

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
def dgmat(mat_float64):
    yield DenseCoancestryMatrix(mat_float64)

def test_is_DenseCoancestryMatrix(dgmat):
    assert is_DenseCoancestryMatrix(dgmat)

def test_ndarray_add_int(dgmat, mat_float64):
    generic_test_ndarray_add(dgmat, 2, mat_float64, 2)

def test_ndarray_add_ndarray(dgmat, mat_float64):
    rmat = numpy.random.binomial(1, 0.5, mat_float64.shape)
    generic_test_ndarray_add(dgmat, rmat, mat_float64, rmat)

def test_ndarray_sub_int(dgmat, mat_float64):
    generic_test_ndarray_sub(dgmat, 2, mat_float64, 2)

def test_ndarray_sub_ndarray(dgmat, mat_float64):
    rmat = numpy.random.binomial(1, 0.5, mat_float64.shape)
    generic_test_ndarray_sub(dgmat, rmat, mat_float64, rmat)
