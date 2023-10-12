import os
import pandas
import pytest
import numpy
import h5py
from pybrops.model.vmat.DenseThreeWayDHAdditiveGenicVarianceMatrix import DenseThreeWayDHAdditiveGenicVarianceMatrix

from pybrops.test.assert_python import assert_concrete_property, not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method

from pybrops.test.model.vmat.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################

### shape fixtures
@pytest.fixture
def ntaxa():
    yield 5

@pytest.fixture
def ntrait():
    yield 3

@pytest.fixture
def ngroup():
    yield 5

### data for variance matrix

@pytest.fixture
def mat(ntaxa,ntrait):
    out = numpy.random.random((ntaxa,ntaxa,ntaxa,ntrait))
    yield out

@pytest.fixture
def taxa(ntaxa):
    out = numpy.array(["Taxon"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def taxa_grp(ngroup,ntaxa):
    out = numpy.random.randint(0, ngroup, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def trait(ntrait):
    out = numpy.array(["Trait"+str(i) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def vmat(mat, taxa, taxa_grp, trait):
    out = DenseThreeWayDHAdditiveGenicVarianceMatrix(
        mat = mat,
        taxa = taxa,
        taxa_grp = taxa_grp,
        trait = trait,
    )
    out.sort_trait()
    out.group_taxa()
    yield out

############################## Test class docstring ############################
def test_class_docstring():
    assert_docstring(DenseThreeWayDHAdditiveGenicVarianceMatrix)

########################### Test concrete properties ###########################

### mat
def test_mat_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "mat")

def test_mat_fget(vmat):
    assert isinstance(vmat.mat, numpy.ndarray)

def test_mat_fset(vmat):
    with not_raises(Exception):
        vmat.mat = numpy.empty((0,0,0,0), dtype=float)

def test_mat_fset_TypeError(vmat):
    with pytest.raises(TypeError):
        vmat.mat = object()

def test_mat_fset_ValueError(vmat):
    with pytest.raises(ValueError):
        vmat.mat = numpy.empty((0,0,0), dtype=float)
    with pytest.raises(ValueError):
        vmat.mat = numpy.empty((0,0,0,0,0), dtype=float)

### square_axes
def test_square_axes_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "square_axes")

### trait_axis
def test_trait_axis_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "trait_axis")

### epgc
def test_epgc_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "epgc")

### nrecurrent
def test_nrecurrent_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "nrecurrent")

### recurrent_axis
def test_recurrent_axis_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "recurrent_axis")

### nfemale
def test_nfemale_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "nfemale")

### female_axis
def test_female_axis_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "female_axis")

### nmale
def test_nmale_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "nmale")

### male_axis
def test_male_axis_is_concrete():
    assert_concrete_property(DenseThreeWayDHAdditiveGenicVarianceMatrix, "male_axis")

############################# Test concrete methods ############################

### __init__
def test___init___is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "__init__")

### to_pandas
def test_to_pandas_is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "to_pandas")

def test_to_pandas(vmat):
    # define column names
    recurrent_col     = "recurrent_col"
    recurrent_grp_col = "recurrent_grp_col"
    female_col        = "female_col"
    female_grp_col    = "female_grp_col"
    male_col          = "male_col"
    male_grp_col      = "male_grp_col"
    trait_col         = "trait_col"
    variance_col      = "variance_col"
    # make dataframe
    df = vmat.to_pandas(
        recurrent_col     = recurrent_col,
        recurrent_grp_col = recurrent_grp_col,
        female_col        = female_col,
        female_grp_col    = female_grp_col,
        male_col          = male_col,
        male_grp_col      = male_grp_col,
        trait_col         = trait_col,
        variance_col      = variance_col,
    )
    # tests
    assert isinstance(df, pandas.DataFrame)
    assert len(df) == vmat.mat.size
    assert len(df.columns) == 8
    assert recurrent_col in df.columns
    assert recurrent_grp_col in df.columns
    assert female_col in df.columns
    assert female_grp_col in df.columns
    assert male_col in df.columns
    assert male_grp_col in df.columns
    assert trait_col in df.columns
    assert variance_col in df.columns
    
    # define column names
    recurrent_col     = "recurrent_col"
    recurrent_grp_col = None
    female_col        = "female_col"
    female_grp_col    = None
    male_col          = "male_col"
    male_grp_col      = None
    trait_col         = "trait_col"
    variance_col      = "variance_col"
    # make dataframe
    df = vmat.to_pandas(
        recurrent_col     = recurrent_col,
        recurrent_grp_col = recurrent_grp_col,
        female_col        = female_col,
        female_grp_col    = female_grp_col,
        male_col          = male_col,
        male_grp_col      = male_grp_col,
        trait_col         = trait_col,
        variance_col      = variance_col,
    )
    # tests
    assert isinstance(df, pandas.DataFrame)
    assert len(df) == vmat.mat.size
    assert len(df.columns) == 5
    assert recurrent_col in df.columns
    assert recurrent_grp_col not in df.columns
    assert female_col in df.columns
    assert female_grp_col not in df.columns
    assert male_col in df.columns
    assert male_grp_col not in df.columns
    assert trait_col in df.columns
    assert variance_col in df.columns

### to_csv
def test_to_csv_is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "to_csv")

def test_to_csv(vmat):
    filename = "test_3dh_vmat.csv"
    # define column names
    recurrent_col     = "recurrent_col"
    recurrent_grp_col = "recurrent_grp_col"
    female_col        = "female_col"
    female_grp_col    = "female_grp_col"
    male_col          = "male_col"
    male_grp_col      = "male_grp_col"
    trait_col         = "trait_col"
    variance_col      = "variance_col"
    # write to file
    vmat.to_csv(
        filename          = filename,
        recurrent_col     = recurrent_col,
        recurrent_grp_col = recurrent_grp_col,
        female_col        = female_col,
        female_grp_col    = female_grp_col,
        male_col          = male_col,
        male_grp_col      = male_grp_col,
        trait_col         = trait_col,
        variance_col      = variance_col,
    )
    # tests
    assert os.path.isfile(filename)

### to_hdf5
def test_to_hdf5_is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "to_hdf5")

def test_to_hdf5(vmat):
    filename = "test_3dh_vmat.h5"
    vmat.to_hdf5(filename)
    assert os.path.isfile(filename)
    # open the file
    h5file = h5py.File(filename, "r")
    if vmat.mat is not None:
        assert "mat" in h5file
    if vmat.taxa is not None:
        assert "taxa" in h5file
    if vmat.taxa_grp is not None:
        assert "taxa_grp" in h5file
    if vmat.trait is not None:
        assert "trait" in h5file
    if vmat.taxa_grp_name is not None:
        assert "taxa_grp_name" in h5file
    if vmat.taxa_grp_stix is not None:
        assert "taxa_grp_stix" in h5file
    if vmat.taxa_grp_spix is not None:
        assert "taxa_grp_spix" in h5file
    if vmat.taxa_grp_len is not None:
        assert "taxa_grp_len" in h5file

### from_pandas
def test_from_pandas_is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "from_pandas")

def test_from_pandas(vmat, ntaxa, ntrait):
    # export
    df = vmat.to_pandas()
    
    # tests for recovery
    out = DenseThreeWayDHAdditiveGenicVarianceMatrix.from_pandas(df)
    assert isinstance(out, DenseThreeWayDHAdditiveGenicVarianceMatrix)
    assert out.ntaxa == ntaxa
    assert out.ntrait == ntrait
    assert out.mat.shape == vmat.mat.shape
    assert out.taxa is not None
    assert all(e in out.taxa for e in vmat.taxa)
    assert all(e in vmat.taxa for e in out.taxa)
    assert out.taxa_grp is not None
    assert all(e in out.taxa_grp for e in vmat.taxa_grp)
    assert all(e in vmat.taxa_grp for e in out.taxa_grp)
    assert numpy.all(out.mat == vmat.mat)

    # tests for recovery without group columns
    out = DenseThreeWayDHAdditiveGenicVarianceMatrix.from_pandas(
        df, 
        recurrent_grp_col = None,
        female_grp_col = None, 
        male_grp_col = None,
    )
    assert isinstance(out, DenseThreeWayDHAdditiveGenicVarianceMatrix)
    out.sort_trait()
    out.group_taxa()
    assert out.ntaxa == ntaxa
    assert out.ntrait == ntrait
    assert out.mat.shape == vmat.mat.shape
    assert out.taxa is not None
    assert all(e in out.taxa for e in vmat.taxa)
    assert all(e in vmat.taxa for e in out.taxa)
    assert out.taxa_grp is None
    assert numpy.all(numpy.isclose(out.mat, vmat.mat))

### from_csv
def test_from_csv_is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "from_csv")

def test_from_csv(vmat, ntaxa, ntrait):
    filename = "test_3dh_vmat.csv"
    # define column names
    recurrent_col     = "recurrent_col"
    recurrent_grp_col = "recurrent_grp_col"
    female_col        = "female_col"
    female_grp_col    = "female_grp_col"
    male_col          = "male_col"
    male_grp_col      = "male_grp_col"
    trait_col         = "trait_col"
    variance_col      = "variance_col"
    # write to file
    vmat.to_csv(
        filename          = filename,
        recurrent_col     = recurrent_col,
        recurrent_grp_col = recurrent_grp_col,
        female_col        = female_col,
        female_grp_col    = female_grp_col,
        male_col          = male_col,
        male_grp_col      = male_grp_col,
        trait_col         = trait_col,
        variance_col      = variance_col,
    )
    # read from file
    out = DenseThreeWayDHAdditiveGenicVarianceMatrix.from_csv(
        filename          = filename,
        recurrent_col     = recurrent_col,
        recurrent_grp_col = recurrent_grp_col,
        female_col        = female_col,
        female_grp_col    = female_grp_col,
        male_col          = male_col,
        male_grp_col      = male_grp_col,
        trait_col         = trait_col,
        variance_col      = variance_col,
    )
    # tests
    assert isinstance(out, DenseThreeWayDHAdditiveGenicVarianceMatrix)
    out.sort_trait()
    out.group_taxa()
    assert out.ntaxa == ntaxa
    assert out.ntrait == ntrait
    assert out.mat.shape == vmat.mat.shape
    assert out.taxa is not None
    assert all(e in out.taxa for e in vmat.taxa)
    assert all(e in vmat.taxa for e in out.taxa)
    assert out.taxa_grp is not None
    assert all(e in out.taxa_grp for e in vmat.taxa_grp)
    assert all(e in vmat.taxa_grp for e in out.taxa_grp)
    assert numpy.all(numpy.isclose(out.mat, vmat.mat))

### from_hdf5
def test_from_hdf5_is_concrete():
    assert_concrete_method(DenseThreeWayDHAdditiveGenicVarianceMatrix, "from_hdf5")

def test_from_hdf5(vmat):
    filename = "test_3dh_vmat.h5"
    vmat.to_hdf5(filename)
    tmp = DenseThreeWayDHAdditiveGenicVarianceMatrix.from_hdf5(filename)
    assert numpy.all(tmp.mat == vmat.mat)
    assert numpy.all(tmp.taxa == vmat.taxa)
    assert numpy.all(tmp.taxa_grp == vmat.taxa_grp)
    assert numpy.all(tmp.taxa_grp_name == vmat.taxa_grp_name)
    assert numpy.all(tmp.taxa_grp_stix == vmat.taxa_grp_stix)
    assert numpy.all(tmp.taxa_grp_spix == vmat.taxa_grp_spix)
    assert numpy.all(tmp.taxa_grp_len == vmat.taxa_grp_len)

