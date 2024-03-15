from pathlib import Path
import pytest
import numpy
import os
import h5py
from os.path import isfile
from pybrops.test.assert_python import assert_classmethod_isabstract, assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.cmat.DenseCoancestryMatrix import check_is_DenseCoancestryMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
################# General shape parameters #################

@pytest.fixture
def ntaxa():
    yield 10 # must be divisible by 2

@pytest.fixture
def nvrnt():
    yield 100

######################## Genotypes #########################
@pytest.fixture
def Xmat(ntaxa, nvrnt):
    out = numpy.random.randint(-1, 2, (ntaxa, nvrnt))
    out = out.astype("int8")
    yield out

@pytest.fixture
def Amat(Xmat):
    A = ((1.0/Xmat.shape[1]) * (Xmat.dot(Xmat.T))) + 1.0
    yield A

@pytest.fixture
def taxa(ntaxa):
    out = numpy.array(["Taxa"+str(i).zfill(3) for i in range(ntaxa)], dtype = object)
    yield out

@pytest.fixture
def taxa_grp(ntaxa):
    out = numpy.repeat(numpy.arange(ntaxa // 2), 2)
    yield out

@pytest.fixture
def cmat(Amat, taxa, taxa_grp):
    out = DummyDenseCoancestryMatrix(
        mat = Amat,
        taxa = taxa,
        taxa_grp = taxa_grp,
    )
    out.group_taxa()
    yield out

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseCoancestryMatrix)

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_gmat_is_abstract():
    assert_classmethod_isabstract(DenseCoancestryMatrix, "from_gmat")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "__init__")

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################
def test_mat_fget(cmat, Amat):
    assert numpy.all(cmat == Amat)

def test_mat_fset_TypeError(cmat, Amat):
    with pytest.raises(TypeError):
        cmat.mat = list(Amat.flatten())

def test_mat_fset_ValueError(cmat, Amat):
    with pytest.raises(ValueError):
        cmat.mat = Amat.flatten()

def test_mat_fset(cmat, Amat):
    cmat.mat = Amat
    assert numpy.all(cmat.mat == Amat)

def test_mat_fdel(cmat, Amat):
    with pytest.raises(AttributeError):
        del cmat.mat

################################################################################
############################# Test concrete methods ############################
################################################################################

############################################################
#################### Matrix conversion #####################

### mat_asformat
def test_mat_asformat_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "mat_asformat")

def test_mat_asformat_coancestry(cmat, Amat):
    A = cmat.mat_asformat("coancestry")
    assert numpy.all(A == Amat)

def test_mat_asformat_kinship(cmat, Amat):
    K = cmat.mat_asformat("kinship")
    assert numpy.all(K == (0.5 * Amat))

def test_mat_asformat_TypeError(cmat):
    with pytest.raises(TypeError):
        K = cmat.mat_asformat([])

def test_mat_asformat_ValueError(cmat):
    with pytest.raises(ValueError):
        K = cmat.mat_asformat("unknown")

############################################################
################ Coancestry/kinship Methods ################

### coancestry
def test_coancestry_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "coancestry")

def test_coancestry_2tuple(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.coancestry(i,j) == Amat[i,j]

def test_coancestry_1tuple(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i) == Amat[i])

def test_coancestry_row_slice(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i,slice(None)) == Amat[i,:])

def test_coancestry_col_slice(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(slice(None),i) == Amat[:,i])

def test_coancestry_list_tuple(cmat, Amat):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.coancestry(a,b) == Amat[a,b])

### kinship
def test_kinship_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "kinship")

def test_kinship_2tuple(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.kinship(i,j) == (0.5 * Amat[i,j])

def test_kinship_1tuple(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(i) == (0.5 * Amat[i]))

def test_kinship_row_slice(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(i,slice(None)) == (0.5 * Amat[i,:]))

def test_kinship_col_slice(cmat, Amat):
    n = Amat.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(slice(None),i) == (0.5 * Amat[:,i]))

def test_kinship_list_tuple(cmat, Amat):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.kinship(a,b) == (0.5 * Amat[a,b]))

### is_positive_semidefinite
def test_is_positive_semidefinite_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "is_positive_semidefinite")

def test_is_positive_semidefinite(cmat):
    assert cmat.is_positive_semidefinite()

def test_is_positive_semidefinite_eigvaltol(cmat):
    assert cmat.is_positive_semidefinite(-1.0)

### apply_jitter
def test_apply_jitter_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "apply_jitter")

def test_apply_jitter(cmat):
    assert cmat.apply_jitter()

### max_inbreeding
def test_max_inbreeding_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "max_inbreeding")

def test_max_inbreeding_coancestry(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.max_inbreeding("coancestry")
    assert tmp == cmat.mat_asformat("coancestry").max()

def test_max_inbreeding_kinship(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.max_inbreeding("kinship")
    assert tmp == cmat.mat_asformat("kinship").max()

def test_max_inbreeding_TypeError(cmat):
    with pytest.raises(TypeError):
        cmat.max_inbreeding(format = object())

def test_max_inbreeding_ValueError(cmat):
    with pytest.raises(ValueError):
        cmat.max_inbreeding(format = "unknown")

### min_inbreeding
def test_min_inbreeding_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "min_inbreeding")

def test_min_inbreeding_TypeError(cmat):
    with pytest.raises(TypeError):
        cmat.min_inbreeding(format = object())

def test_min_inbreeding_ValueError(cmat):
    with pytest.raises(ValueError):
        cmat.min_inbreeding(format = "unknown")

### inverse
def test_inverse_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "inverse")

def test_inverse_coancestry(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.inverse("coancestry")
    assert numpy.all(tmp == numpy.linalg.inv(cmat.mat_asformat("coancestry")))

def test_inverse_kinship(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.inverse("kinship")
    assert numpy.all(tmp == numpy.linalg.inv(cmat.mat_asformat("kinship")))

def test_inverse_TypeError(cmat):
    with pytest.raises(TypeError):
        cmat.inverse(format = object())

def test_inverse_ValueError(cmat):
    with pytest.raises(ValueError):
        cmat.inverse(format = "unknown")

############################################################
################ Matrix summary statistics #################

### max
def test_max_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "max")

def test_max_coancestry(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.max(format = "coancestry")
    assert tmp == cmat.mat_asformat("coancestry").max()

def test_max_kinship(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.max(format = "kinship")
    assert tmp == cmat.mat_asformat("kinship").max()

### mean
def test_mean_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "mean")

def test_mean_coancestry(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.mean(format = "coancestry")
    assert tmp == cmat.mat_asformat("coancestry").mean()

def test_mean_kinship(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.mean(format = "kinship")
    assert tmp == cmat.mat_asformat("kinship").mean()

### min
def test_min_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "min")

def test_min_coancestry(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.min(format = "coancestry")
    assert tmp == cmat.mat_asformat("coancestry").min()

def test_min_kinship(cmat):
    tmp = None
    with not_raises(Exception):
        tmp = cmat.min(format = "kinship")
    assert tmp == cmat.mat_asformat("kinship").min()

############################################################
##################### Matrix File I/O ######################

### to_pandas
def test_to_pandas_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "to_pandas")

def test_to_pandas(cmat):
    # None values
    cmat.taxa = None
    cmat.taxa_grp = None
    df = cmat.to_pandas()
    assert isinstance(df, pandas.DataFrame)
    assert len(df) == cmat.ntaxa
    assert len(df.columns) == (cmat.ntaxa + 2)
    assert numpy.all(df["taxa"] == numpy.array([str(i) for i in range(cmat.ntaxa)]))
    assert numpy.all(df["taxa_grp"].isna())

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = None
    df = cmat.to_pandas()
    assert isinstance(df, pandas.DataFrame)
    assert len(df) == cmat.ntaxa
    assert len(df.columns) == (cmat.ntaxa + 2)
    assert numpy.all(df["taxa"] == cmat.taxa)
    assert numpy.all(df["taxa_grp"].isna())

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = numpy.arange(cmat.ntaxa)
    df = cmat.to_pandas()
    assert isinstance(df, pandas.DataFrame)
    assert len(df) == cmat.ntaxa
    assert len(df.columns) == (cmat.ntaxa + 2)
    assert numpy.all(df["taxa"] == cmat.taxa)
    assert numpy.all(df["taxa_grp"] == cmat.taxa_grp)

def test_to_pandas_taxa_col(cmat):
    with not_raises(Exception):
        df = cmat.to_pandas(taxa_col = "taxa_col")

def test_to_pandas_taxa_col_TypeError(cmat):
    with pytest.raises(TypeError):
        df = cmat.to_pandas(taxa_col = None)
    with pytest.raises(TypeError):
        df = cmat.to_pandas(taxa_col = object())

def test_to_pandas_taxa_grp_col(cmat):
    with not_raises(Exception):
        df = cmat.to_pandas(taxa_grp_col = None)
    with not_raises(Exception):
        df = cmat.to_pandas(taxa_grp_col = "taxa_grp_col")

def test_to_pandas_taxa_grp_col_TypeError(cmat):
    with pytest.raises(TypeError):
        df = cmat.to_pandas(taxa_grp_col = object())

def test_to_pandas_taxa(cmat):
    with not_raises(Exception):
        df = cmat.to_pandas(taxa = "all")
    with not_raises(Exception):
        df = cmat.to_pandas(taxa = [0,1,2,3])

def test_to_pandas_taxa_TypeError(cmat):
    with pytest.raises(TypeError):
        df = cmat.to_pandas(taxa = None)
    with pytest.raises(TypeError):
        df = cmat.to_pandas(taxa = object())

def test_to_pandas_taxa_ValueError(cmat):
    with pytest.raises(ValueError):
        df = cmat.to_pandas(taxa = ["a","b","c","d"])

### to_csv
def test_to_csv_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "to_csv")

def test_to_csv(cmat):
    filename = "sample_coancestry_matrix.csv"
    cmat.to_csv(filename)
    assert isfile(filename)

### to_hdf5
def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseCoancestryMatrix, "to_hdf5")

def test_to_hdf5(cmat):
    filename = "sample_coancestry_matrix.h5"
    if isfile(filename):
        os.remove(filename)
    cmat.to_hdf5(filename)
    assert isfile(filename)

def test_to_hdf5_str(cmat):
    fp = "tmp.h5"
    cmat.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_Path(cmat):
    fp = Path("tmp.h5")
    cmat.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_h5py_File(cmat):
    fp = "tmp.h5"
    h5file = h5py.File(fp, "a")
    with not_raises(Exception):
        cmat.to_hdf5(h5file)
    h5file.close()
    assert os.path.exists(fp)
    os.remove(fp)

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

############################################################
##################### Matrix File I/O ######################

### from_pandas
def test_from_pandas_is_concrete():
    assert_classmethod_isconcrete(DenseCoancestryMatrix, "from_pandas")

def test_from_pandas(cmat):
    # None values
    cmat.taxa = None
    cmat.taxa_grp = None
    df = cmat.to_pandas()
    tmp = DummyDenseCoancestryMatrix.from_pandas(df)
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == numpy.array([str(i) for i in range(cmat.ntaxa)]))
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = None
    df = cmat.to_pandas()
    tmp = DummyDenseCoancestryMatrix.from_pandas(df)
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == cmat.taxa)
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = numpy.arange(cmat.ntaxa)
    df = cmat.to_pandas()
    tmp = DummyDenseCoancestryMatrix.from_pandas(df)
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == cmat.taxa)
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

def test_from_pandas_taxa_col(cmat):
    # None values
    cmat.taxa = None
    cmat.taxa_grp = None
    df = cmat.to_pandas(taxa_col = "taxa_col")
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa_col = "taxa_col")
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == numpy.array([str(i) for i in range(cmat.ntaxa)]))
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = None
    df = cmat.to_pandas(taxa_col = "taxa_col")
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa_col = "taxa_col")
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == cmat.taxa)
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = numpy.arange(cmat.ntaxa)
    df = cmat.to_pandas(taxa_col = "taxa_col")
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa_col = "taxa_col")
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == cmat.taxa)
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

def test_from_pandas_taxa_grp_col(cmat):
    # None values
    cmat.taxa = None
    cmat.taxa_grp = None
    df = cmat.to_pandas(taxa_grp_col = "taxa_grp_col")
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa_grp_col = "taxa_grp_col")
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == numpy.array([str(i) for i in range(cmat.ntaxa)]))
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = None
    df = cmat.to_pandas(taxa_grp_col = "taxa_grp_col")
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa_grp_col = "taxa_grp_col")
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == cmat.taxa)
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = numpy.arange(cmat.ntaxa)
    df = cmat.to_pandas(taxa_grp_col = "taxa_grp_col")
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa_grp_col = "taxa_grp_col")
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == cmat.ntaxa
    assert numpy.all(tmp.mat == cmat.mat)
    assert numpy.all(tmp.taxa == cmat.taxa)
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

def test_from_pandas_taxa(cmat):
    # None values
    cmat.taxa = None
    cmat.taxa_grp = None
    df = cmat.to_pandas(taxa = [0,1,2,3])
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa = [0,1,2,3])
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == len(df)
    cmat_mat_subset = cmat.mat[[0,1,2,3],:][:,[0,1,2,3]]
    assert numpy.all(tmp.mat == cmat_mat_subset)
    assert numpy.all(tmp.taxa == numpy.array([str(i) for i in [0,1,2,3]]))
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = None
    df = cmat.to_pandas(taxa = [0,1,2,3])
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa = [0,1,2,3])
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == len(df)
    cmat_mat_subset = cmat.mat[[0,1,2,3],:][:,[0,1,2,3]]
    assert numpy.all(tmp.mat == cmat_mat_subset)
    assert numpy.all(tmp.taxa == cmat.taxa[[0,1,2,3]])
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp)

    # non-None values
    cmat.taxa = numpy.array([str(i) for i in range(cmat.ntaxa)], dtype=object)
    cmat.taxa_grp = numpy.arange(cmat.ntaxa)
    df = cmat.to_pandas(taxa = [0,1,2,3])
    tmp = DummyDenseCoancestryMatrix.from_pandas(df, taxa = [0,1,2,3])
    assert isinstance(tmp, DenseCoancestryMatrix)
    assert tmp.ntaxa == len(df)
    cmat_mat_subset = cmat.mat[[0,1,2,3],:][:,[0,1,2,3]]
    assert numpy.all(tmp.mat == cmat_mat_subset)
    assert numpy.all(tmp.taxa == cmat.taxa[[0,1,2,3]])
    assert numpy.all(tmp.taxa_grp == cmat.taxa_grp[[0,1,2,3]])

### from_csv
def test_from_csv_is_concrete():
    assert_classmethod_isconcrete(DenseCoancestryMatrix, "from_csv")

def test_from_csv(cmat):
    filename = "saved_coancestry_matrix.csv"
    cmat.to_csv(filename)
    tmp = DummyDenseCoancestryMatrix.from_csv(filename)
    assert isinstance(tmp, DenseCoancestryMatrix)

### from_hdf5
def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DenseCoancestryMatrix, "from_hdf5")

def test_from_hdf5_str(cmat):
    fp = "tmp.h5"
    cmat.to_hdf5(fp)
    out = DummyDenseCoancestryMatrix.from_hdf5(fp)
    # general
    assert numpy.all(cmat.mat == out.mat)
    assert cmat.mat_ndim == out.mat_ndim
    assert cmat.mat_shape == out.mat_shape
    # taxa
    assert numpy.all(cmat.taxa == out.taxa)
    assert numpy.all(cmat.taxa_grp == out.taxa_grp)
    assert cmat.ntaxa == out.ntaxa
    assert cmat.taxa_axis == out.taxa_axis
    assert numpy.all(cmat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(cmat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(cmat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(cmat.taxa_grp_len == out.taxa_grp_len)
    os.remove(fp)

def test_from_hdf5_Path(cmat):
    fp = Path("tmp.h5")
    cmat.to_hdf5(fp)
    out = DummyDenseCoancestryMatrix.from_hdf5(fp)
    # general
    assert numpy.all(cmat.mat == out.mat)
    assert cmat.mat_ndim == out.mat_ndim
    assert cmat.mat_shape == out.mat_shape
    # taxa
    assert numpy.all(cmat.taxa == out.taxa)
    assert numpy.all(cmat.taxa_grp == out.taxa_grp)
    assert cmat.ntaxa == out.ntaxa
    assert cmat.taxa_axis == out.taxa_axis
    assert numpy.all(cmat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(cmat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(cmat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(cmat.taxa_grp_len == out.taxa_grp_len)
    os.remove(fp)

def test_from_hdf5_h5py_File(cmat):
    fp = Path("tmp.h5")
    cmat.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DummyDenseCoancestryMatrix.from_hdf5(h5file)
    # general
    assert numpy.all(cmat.mat == out.mat)
    assert cmat.mat_ndim == out.mat_ndim
    assert cmat.mat_shape == out.mat_shape
    # taxa
    assert numpy.all(cmat.taxa == out.taxa)
    assert numpy.all(cmat.taxa_grp == out.taxa_grp)
    assert cmat.ntaxa == out.ntaxa
    assert cmat.taxa_axis == out.taxa_axis
    assert numpy.all(cmat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(cmat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(cmat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(cmat.taxa_grp_len == out.taxa_grp_len)
    h5file.close()
    os.remove(fp)

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_DenseCoancestryMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseCoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseCoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_DenseCoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_DenseCoancestryMatrix(None, "cmat")
