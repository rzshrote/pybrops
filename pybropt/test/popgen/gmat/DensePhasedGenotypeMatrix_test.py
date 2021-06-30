import copy
import numpy
import pytest
import os.path

from pybropt.popgen.gmat import DensePhasedGenotypeMatrix
from pybropt.popgen.gmat import is_DensePhasedGenotypeMatrix

from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction

################################################################################
############################ Sample Test Variables #############################
################################################################################

@pytest.fixture
def dpgmat(shared_datadir):
    data_path = shared_datadir / "sample.vcf"
    yield DensePhasedGenotypeMatrix.from_vcf(data_path)

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

################################################################################
################################# Sample Tests #################################
################################################################################

def test_is_DensePhasedGenotypeMatrix(dpgmat):
    assert is_DensePhasedGenotypeMatrix(dpgmat)

def test_mat_fget(dpgmat, mat_int8):
    assert numpy.all(dpgmat.mat == mat_int8)

def test_vrnt_chrgrp_fget(dpgmat, mat_chrgrp):
    assert numpy.all(dpgmat.vrnt_chrgrp == mat_chrgrp)

def test_vrnt_phypos_fget(dpgmat, mat_phypos):
    assert numpy.all(dpgmat.vrnt_phypos == mat_phypos)

def test_vrnt_taxa_fget(dpgmat, mat_taxa):
    assert numpy.all(dpgmat.taxa == mat_taxa)

def test_ploidy_fget(dpgmat):
    assert dpgmat.ploidy == dpgmat.mat.shape[0]

def test_nphase_fget(dpgmat):
    assert dpgmat.nphase == dpgmat.mat.shape[0]

def test_ntaxa_fget(dpgmat):
    assert dpgmat.ntaxa == dpgmat.mat.shape[1]

def test_nvrnt_fget(dpgmat):
    assert dpgmat.nvrnt == dpgmat.mat.shape[2]

########################################################
######### Matrix element copy-on-manipulation ##########
########################################################
def test_adjoin_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp

    out = dpgmat.adjoin(dpgmat, axis = 1)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 1)
    out_taxa = numpy.append(mat_taxa, mat_taxa, axis = 0)
    out_taxa_grp = numpy.append(mat_taxa_grp, mat_taxa_grp, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_adjoin_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos

    out = dpgmat.adjoin(dpgmat, axis = 2)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 2)
    out_chrgrp = numpy.append(mat_chrgrp, mat_chrgrp, axis = 0)
    out_phypos = numpy.append(mat_phypos, mat_phypos, axis = 0)
    out_genpos = numpy.append(mat_genpos, mat_genpos, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_adjoin_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp

    out = dpgmat.adjoin_taxa(dpgmat)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 1)
    out_taxa = numpy.append(mat_taxa, mat_taxa, axis = 0)
    out_taxa_grp = numpy.append(mat_taxa_grp, mat_taxa_grp, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_adjoin_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos

    out = dpgmat.adjoin_vrnt(dpgmat)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 2)
    out_chrgrp = numpy.append(mat_chrgrp, mat_chrgrp, axis = 0)
    out_phypos = numpy.append(mat_phypos, mat_phypos, axis = 0)
    out_genpos = numpy.append(mat_genpos, mat_genpos, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_delete_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    out = dpgmat.delete(ix, axis = 1)

    out_mat = numpy.delete(mat_int8, ix, axis = 1)
    out_taxa = numpy.delete(mat_taxa, ix, axis = 0)
    out_taxa_grp = numpy.delete(mat_taxa_grp, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_delete_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1,2])

    out = dpgmat.delete(ix, axis = 2)

    out_mat = numpy.delete(mat_int8, ix, axis = 2)
    out_chrgrp = numpy.delete(mat_chrgrp, ix, axis = 0)
    out_phypos = numpy.delete(mat_phypos, ix, axis = 0)
    out_genpos = numpy.delete(mat_genpos, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_delete_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    out = dpgmat.delete_taxa(ix)

    out_mat = numpy.delete(mat_int8, ix, axis = 1)
    out_taxa = numpy.delete(mat_taxa, ix, axis = 0)
    out_taxa_grp = numpy.delete(mat_taxa_grp, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_delete_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1,2])

    out = dpgmat.delete_vrnt(ix)

    out_mat = numpy.delete(mat_int8, ix, axis = 2)
    out_chrgrp = numpy.delete(mat_chrgrp, ix, axis = 0)
    out_phypos = numpy.delete(mat_phypos, ix, axis = 0)
    out_genpos = numpy.delete(mat_genpos, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_insert_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    out = dpgmat.insert(ix, dpgmat, axis = 1)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 1)
    out_taxa = numpy.insert(mat_taxa, ix, mat_taxa, axis = 0)
    out_taxa_grp = numpy.insert(mat_taxa_grp, ix, mat_taxa_grp, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_insert_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1])

    out = dpgmat.insert(ix, dpgmat, axis = 2)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 2)
    out_chrgrp = numpy.insert(mat_chrgrp, ix, mat_chrgrp, axis = 0)
    out_phypos = numpy.insert(mat_phypos, ix, mat_phypos, axis = 0)
    out_genpos = numpy.insert(mat_genpos, ix, mat_genpos, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_insert_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    out = dpgmat.insert_taxa(ix, dpgmat)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 1)
    out_taxa = numpy.insert(mat_taxa, ix, mat_taxa, axis = 0)
    out_taxa_grp = numpy.insert(mat_taxa_grp, ix, mat_taxa_grp, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_insert_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1])

    out = dpgmat.insert_vrnt(ix, dpgmat)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 2)
    out_chrgrp = numpy.insert(mat_chrgrp, ix, mat_chrgrp, axis = 0)
    out_phypos = numpy.insert(mat_phypos, ix, mat_phypos, axis = 0)
    out_genpos = numpy.insert(mat_genpos, ix, mat_genpos, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_select_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1,2])

    out = dpgmat.select(ix, axis = 1)

    out_mat = numpy.take(mat_int8, ix, axis = 1)
    out_taxa = numpy.take(mat_taxa, ix, axis = 0)
    out_taxa_grp = numpy.take(mat_taxa_grp, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_select_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1,2])

    out = dpgmat.select(ix, axis = 2)

    out_mat = numpy.take(mat_int8, ix, axis = 2)
    out_chrgrp = numpy.take(mat_chrgrp, ix, axis = 0)
    out_phypos = numpy.take(mat_phypos, ix, axis = 0)
    out_genpos = numpy.take(mat_genpos, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_select_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1,2])

    out = dpgmat.select_taxa(ix)

    out_mat = numpy.take(mat_int8, ix, axis = 1)
    out_taxa = numpy.take(mat_taxa, ix, axis = 0)
    out_taxa_grp = numpy.take(mat_taxa_grp, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_select_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1,2])

    out = dpgmat.select_vrnt(ix)

    out_mat = numpy.take(mat_int8, ix, axis = 2)
    out_chrgrp = numpy.take(mat_chrgrp, ix, axis = 0)
    out_phypos = numpy.take(mat_phypos, ix, axis = 0)
    out_genpos = numpy.take(mat_genpos, ix, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_concat_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp

    out = dpgmat.concat([dpgmat, dpgmat], axis = 1)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 1)
    out_taxa = numpy.append(mat_taxa, mat_taxa, axis = 0)
    out_taxa_grp = numpy.append(mat_taxa_grp, mat_taxa_grp, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_concat_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos

    out = dpgmat.concat([dpgmat, dpgmat], axis = 2)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 2)
    out_chrgrp = numpy.append(mat_chrgrp, mat_chrgrp, axis = 0)
    out_phypos = numpy.append(mat_phypos, mat_phypos, axis = 0)
    out_genpos = numpy.append(mat_genpos, mat_genpos, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

def test_concat_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp

    out = dpgmat.concat_taxa([dpgmat, dpgmat])

    out_mat = numpy.append(mat_int8, mat_int8, axis = 1)
    out_taxa = numpy.append(mat_taxa, mat_taxa, axis = 0)
    out_taxa_grp = numpy.append(mat_taxa_grp, mat_taxa_grp, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.taxa == out_taxa)
    assert numpy.all(out.taxa_grp == out_taxa_grp)

def test_concat_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos

    out = dpgmat.concat_vrnt([dpgmat, dpgmat])

    out_mat = numpy.append(mat_int8, mat_int8, axis = 2)
    out_chrgrp = numpy.append(mat_chrgrp, mat_chrgrp, axis = 0)
    out_phypos = numpy.append(mat_phypos, mat_phypos, axis = 0)
    out_genpos = numpy.append(mat_genpos, mat_genpos, axis = 0)

    assert numpy.all(out.mat == out_mat)
    assert numpy.all(out.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(out.vrnt_phypos == out_phypos)
    assert numpy.all(out.vrnt_genpos == out_genpos)

########################################################
######### Matrix element in-place-manipulation #########
########################################################
def test_append_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp

    dpgmat.append(dpgmat, axis = 1)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 1)
    out_taxa = numpy.append(mat_taxa, mat_taxa, axis = 0)
    out_taxa_grp = numpy.append(mat_taxa_grp, mat_taxa_grp, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.taxa == out_taxa)
    assert numpy.all(dpgmat.taxa_grp == out_taxa_grp)

def test_append_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos

    dpgmat.append(dpgmat, axis = 2)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 2)
    out_chrgrp = numpy.append(mat_chrgrp, mat_chrgrp, axis = 0)
    out_phypos = numpy.append(mat_phypos, mat_phypos, axis = 0)
    out_genpos = numpy.append(mat_genpos, mat_genpos, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == out_phypos)
    assert numpy.all(dpgmat.vrnt_genpos == out_genpos)

def test_append_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp

    dpgmat.append_taxa(dpgmat)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 1)
    out_taxa = numpy.append(mat_taxa, mat_taxa, axis = 0)
    out_taxa_grp = numpy.append(mat_taxa_grp, mat_taxa_grp, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.taxa == out_taxa)
    assert numpy.all(dpgmat.taxa_grp == out_taxa_grp)

def test_append_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos

    dpgmat.append_vrnt(dpgmat)

    out_mat = numpy.append(mat_int8, mat_int8, axis = 2)
    out_chrgrp = numpy.append(mat_chrgrp, mat_chrgrp, axis = 0)
    out_phypos = numpy.append(mat_phypos, mat_phypos, axis = 0)
    out_genpos = numpy.append(mat_genpos, mat_genpos, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == out_phypos)
    assert numpy.all(dpgmat.vrnt_genpos == out_genpos)

def test_remove_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    dpgmat.remove(ix, axis = 1)

    out_mat = numpy.delete(mat_int8, ix, axis = 1)
    out_taxa = numpy.delete(mat_taxa, ix, axis = 0)
    out_taxa_grp = numpy.delete(mat_taxa_grp, ix, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.taxa == out_taxa)
    assert numpy.all(dpgmat.taxa_grp == out_taxa_grp)

def test_remove_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1,2])

    dpgmat.remove(ix, axis = 2)

    out_mat = numpy.delete(mat_int8, ix, axis = 2)
    out_chrgrp = numpy.delete(mat_chrgrp, ix, axis = 0)
    out_phypos = numpy.delete(mat_phypos, ix, axis = 0)
    out_genpos = numpy.delete(mat_genpos, ix, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == out_phypos)
    assert numpy.all(dpgmat.vrnt_genpos == out_genpos)

def test_remove_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    dpgmat.remove_taxa(ix)

    out_mat = numpy.delete(mat_int8, ix, axis = 1)
    out_taxa = numpy.delete(mat_taxa, ix, axis = 0)
    out_taxa_grp = numpy.delete(mat_taxa_grp, ix, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.taxa == out_taxa)
    assert numpy.all(dpgmat.taxa_grp == out_taxa_grp)

def test_remove_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1,2])

    dpgmat.remove_vrnt(ix)

    out_mat = numpy.delete(mat_int8, ix, axis = 2)
    out_chrgrp = numpy.delete(mat_chrgrp, ix, axis = 0)
    out_phypos = numpy.delete(mat_phypos, ix, axis = 0)
    out_genpos = numpy.delete(mat_genpos, ix, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == out_phypos)
    assert numpy.all(dpgmat.vrnt_genpos == out_genpos)

def test_incorp_axis_1(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    dpgmat.incorp(ix, dpgmat, axis = 1)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 1)
    out_taxa = numpy.insert(mat_taxa, ix, mat_taxa, axis = 0)
    out_taxa_grp = numpy.insert(mat_taxa_grp, ix, mat_taxa_grp, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.taxa == out_taxa)
    assert numpy.all(dpgmat.taxa_grp == out_taxa_grp)

def test_incorp_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1])

    dpgmat.incorp(ix, dpgmat, axis = 2)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 2)
    out_chrgrp = numpy.insert(mat_chrgrp, ix, mat_chrgrp, axis = 0)
    out_phypos = numpy.insert(mat_phypos, ix, mat_phypos, axis = 0)
    out_genpos = numpy.insert(mat_genpos, ix, mat_genpos, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == out_phypos)
    assert numpy.all(dpgmat.vrnt_genpos == out_genpos)

def test_incorp_taxa(dpgmat, mat_int8, mat_taxa):
    mat_taxa_grp = numpy.arange(len(dpgmat.taxa))
    dpgmat.taxa_grp = mat_taxa_grp
    ix = numpy.int0([1])

    dpgmat.incorp_taxa(ix, dpgmat)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 1)
    out_taxa = numpy.insert(mat_taxa, ix, mat_taxa, axis = 0)
    out_taxa_grp = numpy.insert(mat_taxa_grp, ix, mat_taxa_grp, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.taxa == out_taxa)
    assert numpy.all(dpgmat.taxa_grp == out_taxa_grp)

def test_incorp_vrnt(dpgmat, mat_int8, mat_chrgrp, mat_phypos, mat_genpos):
    dpgmat.vrnt_genpos = mat_genpos
    ix = numpy.int0([1])

    dpgmat.incorp_vrnt(ix, dpgmat)

    out_mat = numpy.insert(mat_int8, ix, mat_int8, axis = 2)
    out_chrgrp = numpy.insert(mat_chrgrp, ix, mat_chrgrp, axis = 0)
    out_phypos = numpy.insert(mat_phypos, ix, mat_phypos, axis = 0)
    out_genpos = numpy.insert(mat_genpos, ix, mat_genpos, axis = 0)

    assert numpy.all(dpgmat.mat == out_mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == out_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == out_phypos)
    assert numpy.all(dpgmat.vrnt_genpos == out_genpos)

########################################################
################### Sorting Methods ####################
########################################################
def test_lexsort_axis_0_error(dpgmat):
    with pytest.raises(RuntimeError):
        dpgmat.lexsort(axis = 0)

def test_lexsort_axis_1_error(dpgmat):
    dpgmat.taxa = None
    with pytest.raises(RuntimeError):
        dpgmat.lexsort(axis = 1)

def test_lexsort_axis_2_error(dpgmat):
    dpgmat.vrnt_chrgrp = None
    dpgmat.vrnt_phypos = None
    with pytest.raises(RuntimeError):
        dpgmat.lexsort(axis = 2)

def test_lexsort_axis_1(dpgmat, mat_taxa):
    a = dpgmat.lexsort(axis = 1)
    assert numpy.all(a == numpy.arange(len(mat_taxa)))

def test_lexsort_axis_2(dpgmat, mat_chrgrp):
    a = dpgmat.lexsort(axis = 2)
    assert numpy.all(a == numpy.arange(len(mat_chrgrp)))

def test_reorder_axis_0(dpgmat, mat_int8):
    ix = numpy.array([1,0])
    dpgmat.reorder(ix, axis = 0)
    mat = mat_int8[ix,:,:]
    assert numpy.all(dpgmat.mat == mat)

def test_reorder_axis_1(dpgmat, mat_int8, mat_taxa):
    ix = numpy.arange(dpgmat.ntaxa)
    numpy.random.shuffle(ix)
    dpgmat.reorder(ix, axis = 1)
    mat = mat_int8[:,ix,:]
    taxa = mat_taxa[ix]
    assert numpy.all(dpgmat.mat == mat)
    assert numpy.all(dpgmat.taxa == taxa)

def test_reorder_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos):
    ix = numpy.arange(dpgmat.nvrnt)
    numpy.random.shuffle(ix)
    dpgmat.reorder(ix, axis = 2)
    mat = mat_int8[:,:,ix]
    chrgrp = mat_chrgrp[ix]
    phypos = mat_phypos[ix]
    assert numpy.all(dpgmat.mat == mat)
    assert numpy.all(dpgmat.vrnt_chrgrp == chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == phypos)

def test_lexsort_axis_1(dpgmat, mat_int8, mat_taxa):
    dpgmat.sort(axis = 1)
    assert numpy.all(dpgmat.mat == mat_int8)
    assert numpy.all(dpgmat.taxa == mat_taxa)

def test_lexsort_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos):
    dpgmat.sort(axis = 1)
    assert numpy.all(dpgmat.mat == mat_int8)
    assert numpy.all(dpgmat.vrnt_chrgrp == mat_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == mat_phypos)

def test_group_axis_1(dpgmat, mat_int8, mat_taxa):
    dpgmat.group(axis = 1)
    assert numpy.all(dpgmat.mat == mat_int8)
    assert numpy.all(dpgmat.taxa == mat_taxa)

def test_group_axis_2(dpgmat, mat_int8, mat_chrgrp, mat_phypos):
    dpgmat.group(axis = 2)
    assert numpy.all(dpgmat.mat == mat_int8)
    assert numpy.all(dpgmat.vrnt_chrgrp == mat_chrgrp)
    assert numpy.all(dpgmat.vrnt_phypos == mat_phypos)

def test_is_grouped_axis_0(dpgmat):
    # dpgmat.group(axis = 0)
    assert not dpgmat.is_grouped(axis = 0)

def test_is_grouped_axis_1(dpgmat):
    dpgmat.taxa_grp = numpy.int64([0,0,1])
    dpgmat.group(axis = 1)
    assert dpgmat.is_grouped(axis = 1)

def test_is_grouped_axis_2(dpgmat):
    dpgmat.group(axis = 2)
    assert dpgmat.is_grouped(axis = 2)

def test_interp_genpos(dpgmat, egmap, mat_genpos):
    dpgmat.interp_genpos(egmap)
    assert numpy.all(dpgmat.vrnt_genpos == mat_genpos)

def test_interp_xoprob_ungrouped(dpgmat, egmap, gmapfn):
    with pytest.raises(RuntimeError):
        dpgmat.interp_xoprob(egmap, gmapfn)

def test_interp_xoprob_grouped(dpgmat, egmap, gmapfn, mat_xoprob):
    dpgmat.group()
    dpgmat.interp_xoprob(egmap, gmapfn)
    assert numpy.allclose(dpgmat.vrnt_xoprob, mat_xoprob)

################################################################################
############################# Song Test Variables ##############################
################################################################################

@pytest.fixture
def song_gmap(shared_datadir):
    data_path = shared_datadir / "Song_2016.linear.M.egmap"
    gmap = ExtendedGeneticMap.from_egmap(data_path)
    gmap.group()
    gmap.build_spline(kind = 'linear', fill_value = 'extrapolate')
    yield gmap

@pytest.fixture
def song_gmapfn():
    yield HaldaneMapFunction()

@pytest.fixture
def song_dpgmat(shared_datadir, song_gmap, song_gmapfn):
    data_path = shared_datadir / "Song_2016_phased_chr_1000.vcf"
    mat = DensePhasedGenotypeMatrix.from_vcf(data_path)
    mat.group()
    mat.interp_xoprob(song_gmap, song_gmapfn)
    yield mat

################################################################################
################################# Song Tests ###################################
################################################################################

def test_xoprob(song_gmap, song_dpgmat):
    assert numpy.all(song_dpgmat.vrnt_xoprob >= 0.0)
    assert numpy.all(song_dpgmat.vrnt_xoprob <= 0.5)

def test_from_to_hdf5(shared_datadir, song_dpgmat):
    # write files
    song_dpgmat.to_hdf5(shared_datadir / "Song_2016_phased_chr_1000.hdf5")
    song_dpgmat.to_hdf5(
        shared_datadir / "Song_2016_phased_chr_1000.hdf5",
        "directoryname"
    )

    # assert file was written
    assert os.path.isfile(shared_datadir / "Song_2016_phased_chr_1000.hdf5")

    # read written files
    song1 = DensePhasedGenotypeMatrix.from_hdf5(
        shared_datadir / "Song_2016_phased_chr_1000.hdf5"
    )
    song2 = DensePhasedGenotypeMatrix.from_hdf5(
        shared_datadir / "Song_2016_phased_chr_1000.hdf5",
        "directoryname"
    )

    # assert file was read correctly
    assert numpy.all(song_dpgmat.mat == song1.mat)
    assert numpy.all(song_dpgmat.vrnt_chrgrp == song1.vrnt_chrgrp)
    assert numpy.all(song_dpgmat.vrnt_phypos == song1.vrnt_phypos)
    assert numpy.all(song_dpgmat.taxa == song1.taxa)
    assert numpy.all(song_dpgmat.taxa_grp == song1.taxa_grp)
    assert numpy.all(song_dpgmat.vrnt_name == song1.vrnt_name)
    assert numpy.all(song_dpgmat.vrnt_genpos == song1.vrnt_genpos)
    assert numpy.all(song_dpgmat.vrnt_xoprob == song1.vrnt_xoprob)
    assert numpy.all(song_dpgmat.vrnt_hapgrp == song1.vrnt_hapgrp)
    assert numpy.all(song_dpgmat.vrnt_mask == song1.vrnt_mask)

    assert numpy.all(song_dpgmat.mat == song2.mat)
    assert numpy.all(song_dpgmat.vrnt_chrgrp == song2.vrnt_chrgrp)
    assert numpy.all(song_dpgmat.vrnt_phypos == song2.vrnt_phypos)
    assert numpy.all(song_dpgmat.taxa == song2.taxa)
    assert numpy.all(song_dpgmat.taxa_grp == song2.taxa_grp)
    assert numpy.all(song_dpgmat.vrnt_name == song2.vrnt_name)
    assert numpy.all(song_dpgmat.vrnt_genpos == song2.vrnt_genpos)
    assert numpy.all(song_dpgmat.vrnt_xoprob == song2.vrnt_xoprob)
    assert numpy.all(song_dpgmat.vrnt_hapgrp == song2.vrnt_hapgrp)
    assert numpy.all(song_dpgmat.vrnt_mask == song2.vrnt_mask)
