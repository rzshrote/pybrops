import numpy
import pytest

from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeMatrix

@pytest.fixture
def mat_int8_unphased():
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
def mat_int8_phased():
    a = numpy.int8([
        [[0, 1, 1, 0, 0, 1, 0, 0],
         [1, 1, 0, 1, 1, 1, 1, 0],
         [1, 0, 0, 0, 0, 1, 1, 0],
         [1, 0, 0, 0, 1, 1, 1, 0],
         [0, 1, 1, 0, 1, 0, 1, 0],
         [0, 0, 1, 1, 0, 1, 0, 1],
         [0, 0, 0, 0, 0, 1, 0, 1],
         [1, 1, 1, 0, 1, 0, 1, 0]],
        [[1, 0, 0, 0, 1, 1, 0, 0],
         [1, 1, 0, 0, 1, 0, 1, 1],
         [0, 0, 0, 1, 0, 0, 1, 0],
         [1, 1, 0, 0, 1, 1, 1, 1],
         [0, 1, 1, 0, 0, 1, 1, 1],
         [1, 1, 0, 1, 0, 1, 1, 0],
         [0, 0, 0, 0, 0, 0, 1, 1],
         [0, 1, 0, 0, 1, 0, 0, 1]]
    ])
    yield a

@pytest.fixture
def dpgmat(mat_int8_phased):
    yield DensePhasedGenotypeMatrix(mat_int8_phased)

def test_is_DensePhasedGenotypeMatrix(dpgmat):
    assert is_DensePhasedGenotypeMatrix(dpgmat)
