import numpy
import pytest

from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeMatrix

@pytest.fixture
def mat_int8():
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
    ]) # (phase, ntaxa, nloci)
    yield a

@pytest.fixture
def dpgmat(mat_int8):
    yield DensePhasedGenotypeMatrix(mat_int8)

def test_is_DensePhasedGenotypeMatrix(dpgmat):
    assert is_DensePhasedGenotypeMatrix(dpgmat)

def test_mat_fset_is_ndarray():
    with pytest.raises(TypeError):
        mat = [[0,1],[1,0]]
        DensePhasedGenotypeMatrix(mat)

def test_mat_fset_ndarray_type():
    with pytest.raises(TypeError):
        mat = numpy.int32([[0,1],[1,0]])
        DensePhasedGenotypeMatrix(mat)

def test_mat_fset_ndarray_ndim():
    with pytest.raises(ValueError):
        mat = numpy.int8([[0,1],[1,0]])
        DensePhasedGenotypeMatrix(mat)

def test_ploidy_fget(dpgmat):
    assert dpgmat.ploidy == dpgmat.mat.shape[0]

def test_nphase_fget(dpgmat):
    assert dpgmat.nphase == dpgmat.mat.shape[0]

def test_ntaxa_fget(dpgmat):
    assert dpgmat.ntaxa == dpgmat.mat.shape[1]

def test_nloci_fget(dpgmat):
    assert dpgmat.nloci == dpgmat.mat.shape[2]
