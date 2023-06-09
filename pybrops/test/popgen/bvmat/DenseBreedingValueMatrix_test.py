import pytest
import numpy
import copy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import check_is_DenseBreedingValueMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_uncentered_float64():
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
def location_float64(mat_uncentered_float64):
    yield mat_uncentered_float64.mean(0)

@pytest.fixture
def scale_float64(mat_uncentered_float64):
    yield mat_uncentered_float64.std(0)

@pytest.fixture
def mat_float64(mat_uncentered_float64, location_float64, scale_float64):
    yield (mat_uncentered_float64 - location_float64) / scale_float64

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
def bvmat(mat_float64, location_float64, scale_float64, taxa_object, taxa_grp_int64, trait_object):
    a = DenseBreedingValueMatrix(
        mat = mat_float64,
        location = location_float64,
        scale = scale_float64,
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
    assert_docstring(DenseBreedingValueMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "__init__")

def test_targmax_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "targmax")

def test_targmin_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "targmin")

def test_tmax_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "tmax")

def test_tmean_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "tmean")

def test_tmin_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "tmin")

def test_trange_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "trange")

def test_tstd_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "tstd")

def test_tvar_is_concrete():
    assert_concrete_method(DenseBreedingValueMatrix, "tvar")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################

################# mat ##################
def test_mat_fget(bvmat, mat_float64):
    assert numpy.all(bvmat.mat == mat_float64)

def test_mat_fset_TypeError(bvmat, mat_float64):
    with pytest.raises(TypeError):
        bvmat.mat = list(mat_float64.flatten())

def test_mat_fset_ValueError(bvmat, mat_float64):
    with pytest.raises(ValueError):
        bvmat.mat = mat_float64.flatten()

def test_mat_fset(bvmat, mat_float64):
    bvmat.mat = mat_float64
    assert numpy.all(bvmat.mat == mat_float64)

def test_mat_fdel(bvmat, mat_float64):
    del bvmat.mat
    with pytest.raises(AttributeError):
        bvmat.mat

############### location ###############
def test_location_fget(bvmat, location_float64):
    assert numpy.all(bvmat.location == location_float64)

def test_location_fset_TypeError(bvmat, location_float64):
    with pytest.raises(TypeError):
        bvmat.location = 5

def test_location_fset_ValueError(bvmat, location_float64):
    l = len(location_float64) // 2
    a = location_float64[0:l]
    with pytest.raises(ValueError):
        bvmat.location = a

def test_location_fset(bvmat, location_float64):
    bvmat.location = location_float64
    assert numpy.all(bvmat.location == location_float64)

def test_location_fdel(bvmat, location_float64):
    del bvmat.location
    with pytest.raises(AttributeError):
        bvmat.location

################ scale #################
def test_scale_fget(bvmat, scale_float64):
    assert numpy.all(bvmat.scale == scale_float64)

def test_scale_fset_TypeError(bvmat, scale_float64):
    with pytest.raises(TypeError):
        bvmat.scale = 5

def test_scale_fset_ValueError(bvmat, scale_float64):
    l = len(scale_float64) // 2
    a = scale_float64[0:l]
    with pytest.raises(ValueError):
        bvmat.scale = a

def test_scale_fset(bvmat, scale_float64):
    bvmat.scale = scale_float64
    assert numpy.all(bvmat.scale == scale_float64)

def test_scale_fdel(bvmat, scale_float64):
    del bvmat.scale
    with pytest.raises(AttributeError):
        bvmat.scale

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_targmax(bvmat):
    a = bvmat.targmax()
    b = numpy.argmax(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_targmin(bvmat):
    a = bvmat.targmin()
    b = numpy.argmin(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_tmax(bvmat):
    a = bvmat.tmax()
    b = numpy.max(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_tmean(bvmat):
    a = bvmat.tmean()
    b = numpy.mean(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_tmin(bvmat):
    a = bvmat.tmin()
    b = numpy.min(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_trange(bvmat):
    a = bvmat.trange()
    b = numpy.ptp(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_tstd(bvmat):
    a = bvmat.tstd()
    b = numpy.std(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

def test_tvar(bvmat):
    a = bvmat.tvar()
    b = numpy.var(bvmat.mat, axis = bvmat.taxa_axis)
    assert numpy.all(a == b)

# TODO: test to_hdf5, from_hdf5

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseBreedingValueMatrix_is_concrete():
    assert_concrete_function(check_is_DenseBreedingValueMatrix)

def test_check_is_DenseBreedingValueMatrix(bvmat):
    with not_raises(TypeError):
        check_is_DenseBreedingValueMatrix(bvmat, "bvmat")
    with pytest.raises(TypeError):
        check_is_DenseBreedingValueMatrix(None, "bvmat")
