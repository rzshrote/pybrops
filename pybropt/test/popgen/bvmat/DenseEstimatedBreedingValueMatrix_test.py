import pytest
import numpy

from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.bvmat import is_DenseEstimatedBreedingValueMatrix

@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [3.3, 9.2, 5.6, 6.8, 8.4],
        [8.7, 3.7, 4.1, 0.4, 1.3],
        [9.0, 4.7, 3.8, 7.9, 3.5]
    ])
    yield a

@pytest.fixture
def debvmat(mat_float64):
    yield DenseEstimatedBreedingValueMatrix(mat_float64)

def test_is_DenseEstimatedBreedingValueMatrix(debvmat):
    assert is_DenseEstimatedBreedingValueMatrix(debvmat)

def test_mat_fget(debvmat, mat_float64):
    assert numpy.all(debvmat.mat == mat_float64)
