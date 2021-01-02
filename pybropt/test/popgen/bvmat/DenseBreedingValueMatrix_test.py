import pytest
import numpy

from pybropt.popgen.bvmat import DenseBreedingValueMatrix
from pybropt.popgen.bvmat import is_DenseBreedingValueMatrix

@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [3.3, 9.2, 5.6, 6.8, 8.4],
        [8.7, 3.7, 4.1, 0.4, 1.3],
        [9.0, 4.7, 3.8, 7.9, 3.5]
    ])
    yield a

@pytest.fixture
def dbvmat(mat_float64):
    yield DenseBreedingValueMatrix(mat_float64)

def test_is_DenseBreedingValueMatrix(dbvmat):
    assert is_DenseBreedingValueMatrix(dbvmat)

def test_mat_fget(dbvmat, mat_float64):
    assert numpy.all(dbvmat.mat == mat_float64)
