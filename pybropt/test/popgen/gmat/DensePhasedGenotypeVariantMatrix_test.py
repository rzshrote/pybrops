import numpy
import pytest

from pybropt.popgen.gmat import DensePhasedGenotypeVariantMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeVariantMatrix

from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction

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
    a = numpy.int64([1, 1, 2, 2, 2, 2, 2, 3])
    yield a

@pytest.fixture
def mat_phypos():
    a = numpy.int64([111, 112, 14370, 17330, 18106, 19302, 19352, 100])
    yield a

@pytest.fixture
def mat_genpos():
    a = numpy.array([0.5275, 0.53, 0.60925, 0.68325, 0.70265, 0.73255, 0.7338, 0.5])
    yield a

@pytest.fixture
def mat_taxa():
    a = numpy.object_(["NA00001", "NA00002", "NA00003"])
    yield a

@pytest.fixture
def mat_xoprob():
    a = numpy.array([0.5, 0.00249376, 0.5, 0.06878444, 0.01902846, 0.02902355, 0.00124844, 0.5])
    yield a

@pytest.fixture
def egmap(shared_datadir):
    e = ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")
    e.group()
    e.build_spline(kind = 'linear', fill_value = 'extrapolate')
    yield e

@pytest.fixture
def gmapfn():
    return HaldaneMapFunction()

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

def test_interp_genpos(dpgvmat, egmap, mat_genpos):
    dpgvmat.interp_genpos(egmap)
    assert numpy.all(dpgvmat.vrnt_genpos == mat_genpos)

def test_interp_xoprob_ungrouped(dpgvmat, egmap, gmapfn):
    with pytest.raises(RuntimeError):
        dpgvmat.interp_xoprob(egmap, gmapfn)

def test_interp_xoprob_grouped(dpgvmat, egmap, gmapfn, mat_xoprob):
    dpgvmat.group()
    dpgvmat.interp_xoprob(egmap, gmapfn)
    assert numpy.allclose(dpgvmat.vrnt_xoprob, mat_xoprob)
