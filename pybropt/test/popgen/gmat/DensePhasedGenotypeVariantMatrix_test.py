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
