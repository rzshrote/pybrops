import pytest
import numpy
import os.path

from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix
from pybropt.popgen.bvmat import is_DenseEstimatedBreedingValueMatrix

################################################################################
############################ Sample Test Variables #############################
################################################################################

@pytest.fixture
def raw_float64():
    numpy.random.seed(453)
    raw = numpy.random.normal(size = (4,5,3))
    yield raw

@pytest.fixture
def mat_float64(raw_float64):
    mat = raw_float64.mean(axis = 0)
    yield mat

@pytest.fixture
def taxa_object_():
    taxa = numpy.object_(["A", "B", "C", "D", "E"])
    yield taxa

@pytest.fixture
def taxa_grp_int64():
    taxa_grp = numpy.int64([1,1,1,2,2])
    yield taxa_grp

@pytest.fixture
def trait_object_():
    trait = numpy.object_(["yield", "oil", "protein"])
    yield trait

@pytest.fixture
def debvmat(mat_float64, raw_float64, taxa_object_, taxa_grp_int64, trait_object_):
    yield DenseEstimatedBreedingValueMatrix(
        mat = mat_float64,
        raw = raw_float64,
        taxa = taxa_object_,
        taxa_grp = taxa_grp_int64,
        trait = trait_object_
    )

################################################################################
################################# Sample Tests #################################
################################################################################

def test_is_DenseEstimatedBreedingValueMatrix(debvmat):
    assert is_DenseEstimatedBreedingValueMatrix(debvmat)

def test_mat_fget(debvmat, mat_float64):
    assert numpy.all(debvmat.mat == mat_float64)

def test_raw_fget(debvmat, raw_float64):
    assert numpy.all(debvmat.raw == raw_float64)

def test_taxa_fget(debvmat, taxa_object_):
    assert numpy.all(debvmat.taxa == taxa_object_)

def test_taxa_grp_fget(debvmat, taxa_grp_int64):
    assert numpy.all(debvmat.taxa_grp == taxa_grp_int64)

def test_trait_fget(debvmat, trait_object_):
    assert numpy.all(debvmat.trait == trait_object_)

########################################################
################### Matrix File I/O ####################
########################################################
def test_from_to_hdf5(shared_datadir, debvmat):
    # write files
    debvmat.to_hdf5(shared_datadir / "test_debvmat.hdf5")
    debvmat.to_hdf5(
        shared_datadir / "test_debvmat.hdf5",
        "directoryname"
    )

    # assert file was written
    assert os.path.isfile(shared_datadir / "test_debvmat.hdf5")

    # read written files
    debvmat1 = DenseEstimatedBreedingValueMatrix.from_hdf5(
        shared_datadir / "test_debvmat.hdf5"
    )
    debvmat2 = DenseEstimatedBreedingValueMatrix.from_hdf5(
        shared_datadir / "test_debvmat.hdf5",
        "directoryname"
    )

    # assert file was read correctly
    assert numpy.all(debvmat.mat == debvmat1.mat)
    assert numpy.all(debvmat.raw == debvmat1.raw)
    assert numpy.all(debvmat.taxa == debvmat1.taxa)
    assert numpy.all(debvmat.taxa_grp == debvmat1.taxa_grp)
    assert numpy.all(debvmat.trait == debvmat1.trait)

    assert numpy.all(debvmat.mat == debvmat2.mat)
    assert numpy.all(debvmat.raw == debvmat2.raw)
    assert numpy.all(debvmat.taxa == debvmat2.taxa)
    assert numpy.all(debvmat.taxa_grp == debvmat2.taxa_grp)
    assert numpy.all(debvmat.trait == debvmat2.trait)
