import pytest
import numpy
import copy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.popgen.bvmat.DenseGenomicEstimatedBreedingValueMatrix import DenseGenomicEstimatedBreedingValueMatrix
from pybrops.popgen.bvmat.DenseGenomicEstimatedBreedingValueMatrix import check_is_DenseGenomicEstimatedBreedingValueMatrix

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
def mat(mat_float64, location_float64, scale_float64, taxa_object, taxa_grp_int64, trait_object):
    a = DenseGenomicEstimatedBreedingValueMatrix(
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
    assert_class_documentation(DenseGenomicEstimatedBreedingValueMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(DenseGenomicEstimatedBreedingValueMatrix, "__init__")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
###################### Test concrete method functionality ######################
################################################################################

# TODO: test to_hdf5, from_hdf5

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseGenomicEstimatedBreedingValueMatrix)

def test_check_is_DenseGenomicEstimatedBreedingValueMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseGenomicEstimatedBreedingValueMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseGenomicEstimatedBreedingValueMatrix(None, "mat")
