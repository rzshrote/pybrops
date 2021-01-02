import pytest
import numpy

from pybropt.popgen.bvmat import DensePhenotypicEstimatedBreedingValueMatrix
from pybropt.popgen.bvmat import is_DensePhenotypicEstimatedBreedingValueMatrix

@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [3.3, 9.2, 5.6, 6.8, 8.4],
        [8.7, 3.7, 4.1, 0.4, 1.3],
        [9.0, 4.7, 3.8, 7.9, 3.5]
    ])
    yield a

@pytest.fixture
def dpebvmat(mat_float64):
    yield DensePhenotypicEstimatedBreedingValueMatrix(mat_float64)

def test_is_DensePhenotypicEstimatedBreedingValueMatrix(dpebvmat):
    assert is_DensePhenotypicEstimatedBreedingValueMatrix(dpebvmat)

def test_mat_fget(dpebvmat, mat_float64):
    assert numpy.all(dpebvmat.mat == mat_float64)
