import pytest
import numpy
import os
from os.path import isfile
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.cmat.DenseCoancestryMatrix import DenseCoancestryMatrix
from pybrops.popgen.cmat.DenseCoancestryMatrix import check_is_DenseCoancestryMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def X_mat_int8():
    X = numpy.array(
       [[-1,  1,  0,  0, -1, -1,  0, -1,  0, -1, -1,  1, -1,  0,  1, -1],
        [-1,  1, -1, -1,  0,  1,  0,  0,  0, -1,  0,  1,  0,  0,  0,  1],
        [ 1,  1,  1,  0,  0,  1,  1, -1, -1,  0,  0,  1, -1, -1, -1, -1],
        [-1,  1,  1, -1,  1,  1,  1,  0,  1,  1,  0,  1,  1, -1, -1,  1],
        [ 1,  1,  1, -1,  0,  0,  0, -1, -1,  1, -1,  0,  1,  0, -1, -1],
        [-1, -1,  1,  0, -1,  0,  1, -1, -1,  0, -1,  1,  0, -1, -1,  1],
        [ 0,  0,  0,  1, -1,  0,  1, -1,  1,  0, -1,  1, -1,  0, -1,  0],
        [-1,  1,  0, -1,  0,  1,  1,  0,  0, -1, -1,  1, -1, -1,  0, -1]], 
        dtype = "int8"
    )
    return X

@pytest.fixture
def A_mat_float64(X_mat_int8):
    A = ((1.0/X_mat_int8.shape[1]) * (X_mat_int8.dot(X_mat_int8.T))) + 1.0
    yield A

@pytest.fixture
def cmat(A_mat_float64):
    yield DummyDenseCoancestryMatrix(A_mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(DenseCoancestryMatrix)

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_gmat_is_abstract():
    assert_abstract_method(DenseCoancestryMatrix, "from_gmat")

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################
def test_mat_fget(cmat, A_mat_float64):
    assert numpy.all(cmat == A_mat_float64)

def test_mat_fset_TypeError(cmat, A_mat_float64):
    with pytest.raises(TypeError):
        cmat.mat = list(A_mat_float64.flatten())

def test_mat_fset_ValueError(cmat, A_mat_float64):
    with pytest.raises(ValueError):
        cmat.mat = A_mat_float64.flatten()

def test_mat_fset(cmat, A_mat_float64):
    cmat.mat = A_mat_float64
    assert numpy.all(cmat.mat == A_mat_float64)

def test_mat_fdel(cmat, A_mat_float64):
    with pytest.raises(AttributeError):
        del cmat.mat

################################################################################
############################# Test concrete methods ############################
################################################################################

### __init__
def test_init_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "__init__")

### mat_asformat
def test_mat_asformat_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "mat_asformat")

def test_mat_asformat_coancestry(cmat, A_mat_float64):
    A = cmat.mat_asformat("coancestry")
    assert numpy.all(A == A_mat_float64)

def test_mat_asformat_kinship(cmat, A_mat_float64):
    K = cmat.mat_asformat("kinship")
    assert numpy.all(K == (0.5 * A_mat_float64))

def test_mat_asformat_TypeError(cmat):
    with pytest.raises(TypeError):
        K = cmat.mat_asformat([])

def test_mat_asformat_ValueError(cmat):
    with pytest.raises(ValueError):
        K = cmat.mat_asformat("unknown")

### coancestry
def test_coancestry_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "coancestry")

def test_coancestry_2tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.coancestry(i,j) == A_mat_float64[i,j]

def test_coancestry_1tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i) == A_mat_float64[i])

def test_coancestry_row_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(i,slice(None)) == A_mat_float64[i,:])

def test_coancestry_col_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.coancestry(slice(None),i) == A_mat_float64[:,i])

def test_coancestry_list_tuple(cmat, A_mat_float64):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.coancestry(a,b) == A_mat_float64[a,b])

### kinship
def test_kinship_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "kinship")

def test_kinship_2tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        for j in range(n):
            assert cmat.kinship(i,j) == (0.5 * A_mat_float64[i,j])

def test_kinship_1tuple(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(i) == (0.5 * A_mat_float64[i]))

def test_kinship_row_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(i,slice(None)) == (0.5 * A_mat_float64[i,:]))

def test_kinship_col_slice(cmat, A_mat_float64):
    n = A_mat_float64.shape[0]
    for i in range(n):
        assert numpy.all(cmat.kinship(slice(None),i) == (0.5 * A_mat_float64[:,i]))

def test_kinship_list_tuple(cmat, A_mat_float64):
    a = [2,3,5]
    b = [1,4,6]
    assert numpy.all(cmat.kinship(a,b) == (0.5 * A_mat_float64[a,b]))

### is_positive_semidefinite
def test_is_positive_semidefinite_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "is_positive_semidefinite")

def test_is_positive_semidefinite(cmat):
    assert cmat.is_positive_semidefinite()

def test_is_positive_semidefinite_eigvaltol(cmat):
    assert cmat.is_positive_semidefinite(-1.0)

### apply_jitter
def test_apply_jitter_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "apply_jitter")

def test_apply_jitter(cmat):
    assert cmat.apply_jitter()

### max_inbreeding
def test_max_inbreeding_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "max_inbreeding")

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
    assert_concrete_method(DenseCoancestryMatrix, "min_inbreeding")

def test_min_inbreeding_TypeError(cmat):
    with pytest.raises(TypeError):
        cmat.min_inbreeding(format = object())

def test_min_inbreeding_ValueError(cmat):
    with pytest.raises(ValueError):
        cmat.min_inbreeding(format = "unknown")

### inverse
def test_inverse_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "inverse")

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

### max
def test_max_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "max")

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
    assert_concrete_method(DenseCoancestryMatrix, "mean")

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
    assert_concrete_method(DenseCoancestryMatrix, "min")

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

### to_pandas
def test_to_pandas_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "to_pandas")

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
    # since taxa values are None
    with pytest.raises(TypeError):
        df = cmat.to_pandas(taxa = ["a","b","c","d"])

def test_to_pandas_taxa_ValueError(cmat):
    cmat.taxa = numpy.array(["Taxon"+str(i) for i in range(cmat.ntaxa)], dtype = object)
    with pytest.raises(ValueError):
        df = cmat.to_pandas(taxa = ["a","b","c","d"])

### to_csv
def test_to_csv_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "to_csv")

def test_to_csv(cmat):
    filename = "sample_coancestry_matrix.csv"
    cmat.to_csv(filename)
    assert isfile(filename)

### to_hdf5
def test_to_hdf5_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "to_hdf5")

def test_to_hdf5(cmat):
    filename = "sample_coancestry_matrix.h5"
    if isfile(filename):
        os.remove(filename)
    cmat.to_hdf5(filename)
    assert isfile(filename)

### from_pandas
def test_from_pandas_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "from_pandas")

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
    assert_concrete_method(DenseCoancestryMatrix, "from_csv")

def test_from_csv(cmat):
    filename = "saved_coancestry_matrix.csv"
    cmat.to_csv(filename)
    tmp = DummyDenseCoancestryMatrix.from_csv(filename)
    assert isinstance(tmp, DenseCoancestryMatrix)

### from_hdf5
def test_from_hdf5_is_concrete():
    assert_concrete_method(DenseCoancestryMatrix, "from_hdf5")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_DenseCoancestryMatrix_is_concrete():
    assert_concrete_function(check_is_DenseCoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseCoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_DenseCoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_DenseCoancestryMatrix(None, "cmat")
