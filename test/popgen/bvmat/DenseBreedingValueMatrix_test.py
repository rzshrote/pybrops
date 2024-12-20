import os
from pathlib import Path
import pandas
import pytest
import numpy
from os.path import isfile
import h5py
from pybrops.test.assert_python import assert_classmethod_isconcrete, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import check_is_DenseBreedingValueMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_uncentered():
    yield numpy.float64([
        [5.9, 5.8, 7. ],
        [5.3, 8.3, 5. ],
        [7.8, 6.4, 7. ],
        [4.8, 7.6, 7.2],
        [5.7, 4.5, 4.8],
        [2. , 7.2, 4.9],
        [5.5, 1.9, 6. ],
        [3.1, 3. , 2.4]
    ])

@pytest.fixture
def location(mat_uncentered):
    yield mat_uncentered.mean(0)

@pytest.fixture
def scale(mat_uncentered):
    yield mat_uncentered.std(0)

@pytest.fixture
def mat(mat_uncentered, location, scale):
    yield (mat_uncentered - location) / scale

###################### Taxa fixtures #######################
@pytest.fixture
def taxa_object():
    a = numpy.object_(["A", "B", "C", "D", "E", "F", "H", "I"])
    yield a

@pytest.fixture
def taxa_grp_int64():
    a = numpy.int64([1,1,2,2,3,3,4,4])
    yield a

@pytest.fixture
def taxa_grp_name_int64():
    a = numpy.int64([1,2,3,4])
    yield a

@pytest.fixture
def taxa_grp_stix_int64():
    a = numpy.int64([0,2,4,6])
    yield a

@pytest.fixture
def taxa_grp_spix_int64():
    a = numpy.int64([2,4,6,8])
    yield a

@pytest.fixture
def taxa_grp_len_int64():
    a = numpy.int64([2,2,2,2])
    yield a

@pytest.fixture
def taxa_lexsort_indices(taxa_object, taxa_grp_int64):
    a = numpy.lexsort((taxa_object, taxa_grp_int64))
    yield a

###################### Trait fixtures ######################
@pytest.fixture
def trait_object():
    a = numpy.object_(["yield", "protein", "oil"])
    yield a

############################################################
@pytest.fixture
def bvmat(mat, location, scale, taxa_object, taxa_grp_int64, trait_object):
    a = DenseBreedingValueMatrix(
        mat = mat,
        location = location,
        scale = scale,
        taxa = taxa_object,
        taxa_grp = taxa_grp_int64,
        trait = trait_object
    )
    a.group_taxa()
    yield a

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseBreedingValueMatrix)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__

def test___init___is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "__init__")

### __copy__
    
def test___copy___is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "__copy__")

### __deepcopy__

def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "__deepcopy__")

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################

### mat
def test_mat_fget(bvmat, mat):
    assert numpy.all(bvmat.mat == mat)

def test_mat_fset(bvmat, mat):
    with not_raises(Exception):
        bvmat.mat = mat
    bvmat.mat = mat
    assert numpy.all(bvmat.mat == mat)

def test_mat_fset_TypeError(bvmat, mat):
    with pytest.raises(TypeError):
        bvmat.mat = object()
    with pytest.raises(TypeError):
        bvmat.mat = list(mat.flatten())

def test_mat_fset_ValueError(bvmat, mat):
    with pytest.raises(ValueError):
        bvmat.mat = mat.flatten()

def test_mat_fdel(bvmat, mat):
    with pytest.raises(AttributeError):
        del bvmat.mat

### location
def test_location_fget(bvmat, location):
    assert numpy.all(bvmat.location == location)

def test_location_fset(bvmat, location):
    with not_raises(Exception):
        bvmat.location = int(5)
    with not_raises(Exception):
        bvmat.location = float(5)
    with not_raises(Exception):
        bvmat.location = location
    bvmat.location = location
    assert numpy.all(bvmat.location == location)

def test_location_fset_TypeError(bvmat, location):
    with pytest.raises(TypeError):
        bvmat.location = object()

def test_location_fset_ValueError(bvmat, location):
    l = len(location) // 2
    a = location[0:l]
    with pytest.raises(ValueError):
        bvmat.location = a

def test_location_fdel(bvmat):
    with pytest.raises(AttributeError):
        del bvmat.location

### scale
def test_scale_fget(bvmat, scale):
    assert numpy.all(bvmat.scale == scale)

def test_scale_fset(bvmat, scale):
    with not_raises(Exception):
        bvmat.scale = int(2)
    with not_raises(Exception):
        bvmat.scale = float(2)
    bvmat.scale = scale
    assert numpy.all(bvmat.scale == scale)

def test_scale_fset_TypeError(bvmat, scale):
    with pytest.raises(TypeError):
        bvmat.scale = object()

def test_scale_fset_ValueError(bvmat, scale):
    with pytest.raises(ValueError):
        bvmat.scale = int(-1)
    with pytest.raises(ValueError):
        bvmat.scale = float(-1)
    with pytest.raises(ValueError):
        l = len(scale) // 2
        a = scale[0:l]
        bvmat.scale = a

def test_scale_fdel(bvmat, scale):
    with pytest.raises(AttributeError):
        del bvmat.scale

################################################################################
############################# Test concrete methods ############################
################################################################################

### targmax

def test_targmax_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "targmax")

def test_targmax(bvmat):
    a = bvmat.targmax()
    b = numpy.argmax(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### targmin

def test_targmin_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "targmin")

def test_targmin(bvmat):
    a = bvmat.targmin()
    b = numpy.argmin(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### tmax

def test_tmax_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "tmax")

def test_tmax(bvmat):
    a = bvmat.tmax()
    b = numpy.max(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### tmean

def test_tmean_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "tmean")

def test_tmean(bvmat):
    a = bvmat.tmean()
    b = numpy.mean(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### tmin

def test_tmin_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "tmin")

def test_tmin(bvmat):
    a = bvmat.tmin()
    b = numpy.min(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### trange

def test_trange_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "trange")

def test_trange(bvmat):
    a = bvmat.trange()
    b = numpy.ptp(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### tstd

def test_tstd_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "tstd")

def test_tstd(bvmat):
    a = bvmat.tstd()
    b = numpy.std(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

### tvar

def test_tvar_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "tvar")

def test_tvar(bvmat):
    a = bvmat.tvar()
    b = numpy.var(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

############################################################
##################### Matrix File I/O ######################

### to_pandas
def test_to_pandas(bvmat, trait_object):
    df = bvmat.to_pandas()
    assert "taxa" in df.columns
    assert "taxa_grp" in df.columns
    for t in trait_object:
        assert t in df.columns

def test_to_pandas_col_rename(bvmat):
    taxa_col = "taxa_col"
    taxa_grp_col = "taxa_grp_col"
    trait_cols = ["t1","t2","t3"]
    df = bvmat.to_pandas(
        taxa_col = taxa_col,
        taxa_grp_col = taxa_grp_col,
        trait_cols = trait_cols
    )
    assert taxa_col in df.columns
    assert taxa_grp_col in df.columns
    for t in trait_cols:
        assert t in df.columns

def test_to_pandas_unscale(bvmat):
    df1 = bvmat.to_pandas(unscale = False)
    df2 = bvmat.to_pandas(unscale = True)
    assert numpy.any(df1 != df2)

def test_to_pandas_TypeError(bvmat):
    with pytest.raises(TypeError):
        df = bvmat.to_pandas(
            taxa_col = object(),
            taxa_grp_col = "taxa_grp",
            trait_cols = ["Trait1","Trait2","Trait3"]
        )
    with pytest.raises(TypeError):
        df = bvmat.to_pandas(
            taxa_col = "taxa",
            taxa_grp_col = object(),
            trait_cols = ["Trait1","Trait2","Trait3"]
        )
    with pytest.raises(TypeError):
        df = bvmat.to_pandas(
            taxa_col = "taxa",
            taxa_grp_col = "taxa_grp",
            trait_cols = object()
        )

def test_to_pandas_ValueError(bvmat):
    taxa_col = "taxa_col"
    taxa_grp_col = "taxa_grp_col"
    trait_cols = ["t"+str(i) for i in range(bvmat.ntrait-1)]
    with pytest.raises(ValueError):
        df = bvmat.to_pandas(
            taxa_col = taxa_col,
            taxa_grp_col = taxa_grp_col,
            trait_cols = trait_cols
        )

### to_csv
def test_to_csv(bvmat):
    filename = "sample_breeding_values.csv"
    bvmat.to_csv(filename)
    assert isfile(filename)

### to_hdf5

def test_to_hdf5_is_concrete():
    assert_method_isconcrete(DenseBreedingValueMatrix, "to_hdf5")

def test_to_hdf5_str(bvmat):
    fp = "tmp.h5"
    bvmat.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_Path(bvmat):
    fp = Path("tmp.h5")
    bvmat.to_hdf5(fp)
    assert os.path.exists(fp)
    os.remove(fp)

def test_to_hdf5_h5py_File(bvmat):
    fp = "tmp.h5"
    h5file = h5py.File(fp, "a")
    with not_raises(Exception):
        bvmat.to_hdf5(h5file)
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
    assert_classmethod_isconcrete(DenseBreedingValueMatrix, "from_pandas")

def test_from_pandas():
    ntaxa = 100
    taxa = numpy.array(["Taxon"+str(i) for i in range(ntaxa)])
    taxa_grp = numpy.random.randint(0,10,ntaxa)
    trait1 = numpy.random.random(ntaxa)
    trait2 = numpy.random.random(ntaxa)
    trait3 = numpy.random.random(ntaxa)
    mat = numpy.stack([trait1,trait2,trait3], axis = 1)
    df = pandas.DataFrame({
        "taxa": taxa,
        "taxa_grp": taxa_grp,
        "trait1": trait1,
        "trait2": trait2,
        "trait3": trait3,
    })
    bvmat = DenseBreedingValueMatrix.from_pandas(
        df = df,
        location = 0.0,
        scale = 1.0,
    )
    assert isinstance(bvmat, DenseBreedingValueMatrix)
    assert numpy.all(bvmat.taxa == taxa)
    assert numpy.all(bvmat.taxa_grp == taxa_grp)
    assert numpy.all(bvmat.mat == mat)

def test_from_pandas_col_rename():
    ntaxa = 100
    taxa = numpy.array(["Taxon"+str(i) for i in range(ntaxa)])
    taxa_grp = numpy.random.randint(0,10,ntaxa)
    trait1 = numpy.random.random(ntaxa)
    trait2 = numpy.random.random(ntaxa)
    trait3 = numpy.random.random(ntaxa)
    mat = numpy.stack([trait1,trait2,trait3], axis = 1)
    df = pandas.DataFrame({
        "taxa_col": taxa,
        "taxa_grp_col": taxa_grp,
        "trait1": trait1,
        "trait2": trait2,
        "trait3": trait3,
    })
    bvmat = DenseBreedingValueMatrix.from_pandas(
        df = df,
        location = 0.0,
        scale = 1.0,
        taxa_col = "taxa_col",
        taxa_grp_col = "taxa_grp_col",
        trait_cols = ["trait1","trait2","trait3"]
    )
    assert isinstance(bvmat, DenseBreedingValueMatrix)
    assert numpy.all(bvmat.taxa == taxa)
    assert numpy.all(bvmat.taxa_grp == taxa_grp)
    assert numpy.all(bvmat.mat == mat)

def test_from_csv():
    filename = "sample_breeding_values.csv"
    out = DenseBreedingValueMatrix.from_csv(filename)
    assert isinstance(out, DenseBreedingValueMatrix)

### from_hdf5

def test_from_hdf5_is_concrete():
    assert_classmethod_isconcrete(DenseBreedingValueMatrix, "from_hdf5")

def test_from_hdf5_str(bvmat):
    fp = "tmp.h5"
    bvmat.to_hdf5(fp)
    out = DenseBreedingValueMatrix.from_hdf5(fp)
    # general
    assert numpy.all(bvmat.mat == out.mat)
    assert bvmat.mat_ndim == out.mat_ndim
    assert bvmat.mat_shape == out.mat_shape
    assert numpy.all(bvmat.location == out.location)
    assert numpy.all(bvmat.scale == out.scale)
    # taxa
    assert numpy.all(bvmat.taxa == out.taxa)
    assert numpy.all(bvmat.taxa_grp == out.taxa_grp)
    assert bvmat.ntaxa == out.ntaxa
    assert bvmat.taxa_axis == out.taxa_axis
    assert numpy.all(bvmat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(bvmat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(bvmat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(bvmat.taxa_grp_len == out.taxa_grp_len)
    # trait
    assert numpy.all(bvmat.trait == out.trait)
    assert bvmat.ntrait == out.ntrait
    assert bvmat.trait_axis == out.trait_axis
    os.remove(fp)

def test_from_hdf5_Path(bvmat):
    fp = Path("tmp.h5")
    bvmat.to_hdf5(fp)
    out = DenseBreedingValueMatrix.from_hdf5(fp)
    # general
    assert numpy.all(bvmat.mat == out.mat)
    assert bvmat.mat_ndim == out.mat_ndim
    assert bvmat.mat_shape == out.mat_shape
    assert numpy.all(bvmat.location == out.location)
    assert numpy.all(bvmat.scale == out.scale)
    # taxa
    assert numpy.all(bvmat.taxa == out.taxa)
    assert numpy.all(bvmat.taxa_grp == out.taxa_grp)
    assert bvmat.ntaxa == out.ntaxa
    assert bvmat.taxa_axis == out.taxa_axis
    assert numpy.all(bvmat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(bvmat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(bvmat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(bvmat.taxa_grp_len == out.taxa_grp_len)
    # trait
    assert numpy.all(bvmat.trait == out.trait)
    assert bvmat.ntrait == out.ntrait
    assert bvmat.trait_axis == out.trait_axis
    os.remove(fp)

def test_from_hdf5_h5py_File(bvmat):
    fp = Path("tmp.h5")
    bvmat.to_hdf5(fp)
    h5file = h5py.File(fp)
    out = DenseBreedingValueMatrix.from_hdf5(h5file)
    # general
    assert numpy.all(bvmat.mat == out.mat)
    assert bvmat.mat_ndim == out.mat_ndim
    assert bvmat.mat_shape == out.mat_shape
    assert numpy.all(bvmat.location == out.location)
    assert numpy.all(bvmat.scale == out.scale)
    # taxa
    assert numpy.all(bvmat.taxa == out.taxa)
    assert numpy.all(bvmat.taxa_grp == out.taxa_grp)
    assert bvmat.ntaxa == out.ntaxa
    assert bvmat.taxa_axis == out.taxa_axis
    assert numpy.all(bvmat.taxa_grp_name == out.taxa_grp_name)
    assert numpy.all(bvmat.taxa_grp_stix == out.taxa_grp_stix)
    assert numpy.all(bvmat.taxa_grp_spix == out.taxa_grp_spix)
    assert numpy.all(bvmat.taxa_grp_len == out.taxa_grp_len)
    # trait
    assert numpy.all(bvmat.trait == out.trait)
    assert bvmat.ntrait == out.ntrait
    assert bvmat.trait_axis == out.trait_axis
    h5file.close()
    os.remove(fp)


################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseBreedingValueMatrix)

def test_check_is_DenseBreedingValueMatrix(bvmat):
    with not_raises(TypeError):
        check_is_DenseBreedingValueMatrix(bvmat, "bvmat")
    with pytest.raises(TypeError):
        check_is_DenseBreedingValueMatrix(None, "bvmat")
