import pytest
import numpy
import copy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.popgen.bvmat import DenseBreedingValueMatrix
from pybropt.popgen.bvmat import is_DenseBreedingValueMatrix
from pybropt.popgen.bvmat import check_is_DenseBreedingValueMatrix
from pybropt.popgen.bvmat import cond_check_is_DenseBreedingValueMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [5.9, 5.8, 7. ],
        [5.3, 8.3, 5. ],
        [7.8, 6.4, 7. ],
        [4.8, 7.6, 7.2],
        [5.7, 4.5, 4.8],
        [2. , 7.2, 4.9],
        [5.5, 1.9, 6. ],
        [3.1, 3. , 2.4]
    ])
    yield a

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
def mat(mat_float64, taxa_object, taxa_grp_int64, trait_object):
    a = DenseBreedingValueMatrix(
        mat = mat_float64,
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
    generic_assert_docstring(DenseBreedingValueMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "__init__")

def test_targmax_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "targmax")

def test_targmin_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "targmin")

def test_tmax_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "tmax")

def test_tmean_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "tmean")

def test_tmin_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "tmin")

def test_trange_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "trange")

def test_tstd_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "tstd")

def test_tvar_is_concrete():
    generic_assert_concrete_method(DenseBreedingValueMatrix, "tvar")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################ General matrix properties #################
def test_mat_fget(mat, mat_float64):
    assert numpy.all(mat == mat_float64)

def test_mat_fset_TypeError(mat, mat_float64):
    with pytest.raises(TypeError):
        mat.mat = list(mat_float64.flatten())

def test_mat_fset_ValueError(mat, mat_float64):
    with pytest.raises(ValueError):
        mat.mat = mat_float64.flatten()

def test_mat_fset(mat, mat_float64):
    mat.mat = mat_float64
    assert numpy.all(mat.mat == mat_float64)

def test_mat_fdel(mat, mat_float64):
    del mat.mat
    with pytest.raises(AttributeError):
        mat.mat

################################################################################
###################### Test concrete method functionality ######################
################################################################################
def test_targmax(mat):
    a = mat.targmax()
    b = numpy.argmax(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_targmin(mat):
    a = mat.targmin()
    b = numpy.argmin(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_tmax(mat):
    a = mat.tmax()
    b = numpy.max(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_tmean(mat):
    a = mat.tmean()
    b = numpy.mean(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_tmin(mat):
    a = mat.tmin()
    b = numpy.min(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_trange(mat):
    a = mat.trange()
    b = numpy.ptp(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_tstd(mat):
    a = mat.tstd()
    b = numpy.std(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

def test_tvar(mat):
    a = mat.tvar()
    b = numpy.var(mat.mat, axis = mat.taxa_axis)
    assert numpy.all(a == b)

# TODO: test to_hdf5, from_hdf5

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseBreedingValueMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseBreedingValueMatrix)

def test_is_DenseBreedingValueMatrix(mat):
    assert is_DenseBreedingValueMatrix(mat)

def test_check_is_DenseBreedingValueMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseBreedingValueMatrix)

def test_check_is_DenseBreedingValueMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseBreedingValueMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseBreedingValueMatrix(None, "mat")

def test_cond_check_is_DenseBreedingValueMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseBreedingValueMatrix)
