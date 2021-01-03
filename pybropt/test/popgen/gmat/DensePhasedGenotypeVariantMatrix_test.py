import numpy
import pytest

from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeVariantMatrix

@pytest.fixture
def dpgvmat(shared_datadir):
    data_path = shared_datadir / "sample.vcf"
    yield DensePhasedGenotypeVariantMatrix.from_vcf(data_path)

@pytest.fixture
def mat_int8():
    a = numpy.int8([
        [[0, 0, 0, 0, 1, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 0, 0],
         [0, 0, 1, 0, 0, 0, 1, 0]],
        [[0, 0, 0, 0, 0, 0, 0, 1],
         [0, 0, 0, 1, 1, 0, 0, 1],
         [1, 1, 1, 0, 0, 0, 1, 1]]
    ])
    yield a

@pytest.fixture
def mat_chrgrp():
    a = numpy.int64([19, 19, 20, 20, 20, 20, 20, 21])
    yield a

@pytest.fixture
def mat_phypos():
    a = numpy.int64([111, 112, 14370, 17330, 1110696, 1230237, 1235237, 10])
    yield a

@pytest.fixture
def mat_taxa():
    a = numpy.string_(["NA00001", "NA00002", "NA00003"])
    yield a

def test_is_DensePhasedGenotypeVariantMatrix(dpgvmat):
    assert is_DensePhasedGenotypeVariantMatrix(dpgvmat)

def test_mat_fget(dpgvmat, mat_int8):
    assert numpy.all(dpgvmat.mat == mat_int8)

def test_vrnt_chrgrp_fget(dpgvmat, mat_chrgrp):
    assert numpy.all(dpgvmat.vrnt_chrgrp == mat_chrgrp)

def test_vrnt_phypos_fget(dpgvmat, mat_phypos):
    assert numpy.all(dpgvmat.vrnt_phypos == mat_phypos)

def test_vrnt_taxa_fget(dpgvmat, mat_taxa):
    assert numpy.all(dpgvmat.taxa == mat_taxa)

def test_axis_index(dpgvmat):
    assert dpgvmat.axis_index(-3) == 0
    assert dpgvmat.axis_index(-2) == 1
    assert dpgvmat.axis_index(-1) == 2
    assert dpgvmat.axis_index(0) == 0
    assert dpgvmat.axis_index(1) == 1
    assert dpgvmat.axis_index(2) == 2

def test_axis_index_error_low(dpgvmat):
    with pytest.raises(IndexError):
        dpgvmat.axis_index(-4)

def test_axis_index_error_high(dpgvmat):
    with pytest.raises(IndexError):
        dpgvmat.axis_index(3)

def test_lexsort_axis_0_error(dpgvmat):
    with pytest.raises(RuntimeError):
        dpgvmat.lexsort(axis = 0)

def test_lexsort_axis_1_error(dpgvmat):
    dpgvmat.taxa = None
    with pytest.raises(RuntimeError):
        dpgvmat.lexsort(axis = 1)

def test_lexsort_axis_2_error(dpgvmat):
    dpgvmat.vrnt_chrgrp = None
    dpgvmat.vrnt_phypos = None
    with pytest.raises(RuntimeError):
        dpgvmat.lexsort(axis = 2)

def test_lexsort_axis_1(dpgvmat, mat_taxa):
    a = dpgvmat.lexsort(axis = 1)
    assert numpy.all(a == numpy.arange(len(mat_taxa)))

def test_lexsort_axis_2(dpgvmat, mat_chrgrp):
    a = dpgvmat.lexsort(axis = 2)
    assert numpy.all(a == numpy.arange(len(mat_chrgrp)))

def test_reorder_axis_0(dpgvmat, mat_int8):
    ix = numpy.array([1,0])
    dpgvmat.reorder(ix, axis = 0)
    mat = mat_int8[ix,:,:]
    assert numpy.all(dpgvmat.mat == mat)

def test_reorder_axis_1(dpgvmat, mat_int8, mat_taxa):
    ix = numpy.arange(dpgvmat.ntaxa)
    numpy.random.shuffle(ix)
    dpgvmat.reorder(ix, axis = 1)
    mat = mat_int8[:,ix,:]
    taxa = mat_taxa[ix]
    assert numpy.all(dpgvmat.mat == mat)
    assert numpy.all(dpgvmat.taxa == taxa)

def test_reorder_axis_2(dpgvmat, mat_int8, mat_chrgrp, mat_phypos):
    ix = numpy.arange(dpgvmat.nloci)
    numpy.random.shuffle(ix)
    dpgvmat.reorder(ix, axis = 2)
    mat = mat_int8[:,:,ix]
    chrgrp = mat_chrgrp[ix]
    phypos = mat_phypos[ix]
    assert numpy.all(dpgvmat.mat == mat)
    assert numpy.all(dpgvmat.vrnt_chrgrp == chrgrp)
    assert numpy.all(dpgvmat.vrnt_phypos == phypos)

def test_lexsort_axis_1(dpgvmat, mat_int8, mat_taxa):
    dpgvmat.sort(axis = 1)
    assert numpy.all(dpgvmat.mat == mat_int8)
    assert numpy.all(dpgvmat.taxa == mat_taxa)

def test_lexsort_axis_2(dpgvmat, mat_int8, mat_chrgrp, mat_phypos):
    dpgvmat.sort(axis = 1)
    assert numpy.all(dpgvmat.mat == mat_int8)
    assert numpy.all(dpgvmat.vrnt_chrgrp == mat_chrgrp)
    assert numpy.all(dpgvmat.vrnt_phypos == mat_phypos)

def test_group_axis_1(dpgvmat, mat_int8, mat_taxa):
    dpgvmat.group(axis = 1)
    assert numpy.all(dpgvmat.mat == mat_int8)
    assert numpy.all(dpgvmat.taxa == mat_taxa)

def test_group_axis_2(dpgvmat, mat_int8, mat_chrgrp, mat_phypos):
    dpgvmat.group(axis = 2)
    assert numpy.all(dpgvmat.mat == mat_int8)
    assert numpy.all(dpgvmat.vrnt_chrgrp == mat_chrgrp)
    assert numpy.all(dpgvmat.vrnt_phypos == mat_phypos)

def test_is_grouped_axis_0(dpgvmat):
    # dpgvmat.group(axis = 0)
    assert not dpgvmat.is_grouped(axis = 0)

def test_is_grouped_axis_1(dpgvmat):
    dpgvmat.taxa_grp = numpy.int64([0,0,1])
    dpgvmat.group(axis = 1)
    assert dpgvmat.is_grouped(axis = 1)

def test_is_grouped_axis_2(dpgvmat):
    dpgvmat.group(axis = 2)
    assert dpgvmat.is_grouped(axis = 2)
