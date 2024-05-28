import numpy
import pytest

from pybrops.model.pmebvmat.DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix import DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix, check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix
from pybrops.test.assert_python import assert_class_isconcrete, assert_classmethod_isconcrete, assert_function_isconcrete, assert_method_isconcrete, assert_module_documentation, assert_module_public_api, not_raises

############################################################
################# General shape parameters #################
############################################################

@pytest.fixture
def ngroup():
    yield 5

@pytest.fixture
def ntaxa(ngroup):
    yield ngroup * 20

@pytest.fixture
def ntrait():
    yield 3

############################################################
####### Progeny Mean Estimated Breeding Value Matrix #######
############################################################

@pytest.fixture
def pmebvmat_mat(ntaxa, ntrait):
    out = numpy.random.random((ntaxa,ntaxa,ntrait))
    yield out

@pytest.fixture
def pmebvmat_location(ntrait):
    out = numpy.repeat(0.0, ntrait)
    yield out

@pytest.fixture
def pmebvmat_scale(ntrait):
    out = numpy.repeat(1.0, ntrait)
    yield out

@pytest.fixture
def pmebvmat_taxa(ntaxa):
    out = numpy.array(["Line"+str(i).zfill(3) for i in range(ntaxa)], dtype=object)
    yield out

@pytest.fixture
def pmebvmat_taxa_grp(ntaxa, ngroup):
    out = numpy.random.randint(1, ngroup+1, ntaxa)
    out.sort()
    yield out

@pytest.fixture
def pmebvmat_trait(ntrait):
    out = numpy.array(["Trait"+str(i).zfill(3) for i in range(ntrait)], dtype=object)
    yield out

@pytest.fixture
def pmebvmat(
        pmebvmat_mat,
        pmebvmat_location,
        pmebvmat_scale,
        pmebvmat_taxa,
        pmebvmat_taxa_grp,
        pmebvmat_trait,
    ):
    out = DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix(
        mat = pmebvmat_mat,
        location = pmebvmat_location,
        scale = pmebvmat_scale,
        taxa = pmebvmat_taxa,
        taxa_grp = pmebvmat_taxa_grp,
        trait = pmebvmat_trait,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.pmebvmat.DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix
    assert_module_documentation(pybrops.model.pmebvmat.DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix)

def test_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.pmebvmat.DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix
    assert_module_public_api(pybrops.model.pmebvmat.DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix)

################################################################################
########################### Test class documentation ###########################
################################################################################

def test_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix_is_concrete():
    assert_class_isconcrete(DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix)

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix, "__init__")

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix
def test_check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix)

def test_check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix(pmebvmat):
    with not_raises(TypeError):
        check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix(pmebvmat, "pmebvmat")
    with pytest.raises(TypeError):
        check_is_DenseTwoWayProgenyMeanEstimatedBreedingValueMatrix(None, "pmebvmat")