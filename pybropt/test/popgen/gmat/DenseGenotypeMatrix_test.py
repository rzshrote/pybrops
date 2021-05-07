import numpy
import pytest

from pybropt.popgen.gmat import DenseGenotypeMatrix
from pybropt.popgen.gmat import is_DenseGenotypeMatrix

from pybropt.test import generic_test_ndarray_add
from pybropt.test import generic_test_ndarray_sub

@pytest.fixture
def mat_int8():
    a = numpy.int8([
        [0, 1, 1, 0, 0, 1, 0, 0],
        [1, 1, 0, 1, 1, 1, 1, 0],
        [1, 0, 0, 0, 0, 1, 1, 0],
        [1, 0, 0, 0, 1, 1, 1, 0],
        [0, 1, 1, 0, 1, 0, 1, 0],
        [0, 0, 1, 1, 0, 1, 0, 1],
        [0, 0, 0, 0, 0, 1, 0, 1],
        [1, 1, 1, 0, 1, 0, 1, 0]
    ])
    yield a

@pytest.fixture
def dgmat(mat_int8):
    yield DenseGenotypeMatrix(mat_int8)

def test_is_DenseGenotypeMatrix(dgmat):
    assert is_DenseGenotypeMatrix(dgmat)

def test_ndarray_add_int(dgmat, mat_int8):
    generic_test_ndarray_add(dgmat, 2, mat_int8, 2)

def test_ndarray_add_ndarray(dgmat, mat_int8):
    rmat = numpy.random.binomial(1, 0.5, mat_int8.shape)
    generic_test_ndarray_add(dgmat, rmat, mat_int8, rmat)

def test_ndarray_sub_int(dgmat, mat_int8):
    generic_test_ndarray_sub(dgmat, 2, mat_int8, 2)

def test_ndarray_sub_ndarray(dgmat, mat_int8):
    rmat = numpy.random.binomial(1, 0.5, mat_int8.shape)
    generic_test_ndarray_sub(dgmat, rmat, mat_int8, rmat)
