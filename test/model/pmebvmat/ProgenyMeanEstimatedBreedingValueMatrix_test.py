import pytest

from pybrops.model.pmebvmat.ProgenyMeanEstimatedBreedingValueMatrix import check_is_ProgenyMeanEstimatedBreedingValueMatrix
from pybrops.test.assert_python import assert_class_isabstract, assert_function_isconcrete, assert_module_documentation, assert_module_public_api, not_raises

from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
####### Progeny Mean Estimated Breeding Value Matrix #######
############################################################

@pytest.fixture
def pmebvmat():
    out = DummyProgenyMeanEstimatedBreedingValueMatrix()
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_ProgenyMeanEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.pmebvmat.ProgenyMeanEstimatedBreedingValueMatrix
    assert_module_documentation(pybrops.model.pmebvmat.ProgenyMeanEstimatedBreedingValueMatrix)

def test_ProgenyMeanEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.pmebvmat.ProgenyMeanEstimatedBreedingValueMatrix
    assert_module_public_api(pybrops.model.pmebvmat.ProgenyMeanEstimatedBreedingValueMatrix)

################################################################################
########################### Test class documentation ###########################
################################################################################

def test_ProgenyMeanEstimatedBreedingValueMatrix_is_abstract():
    assert_class_isabstract(ProgenyMeanEstimatedBreedingValueMatrix)

################################################################################
############################ Test Class Properties #############################
################################################################################

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_ProgenyMeanEstimatedBreedingValueMatrix
def test_check_is_ProgenyMeanEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_ProgenyMeanEstimatedBreedingValueMatrix)

def test_check_is_ProgenyMeanEstimatedBreedingValueMatrix(pmebvmat):
    with not_raises(TypeError):
        check_is_ProgenyMeanEstimatedBreedingValueMatrix(pmebvmat, "pmebvmat")
    with pytest.raises(TypeError):
        check_is_ProgenyMeanEstimatedBreedingValueMatrix(None, "pmebvmat")
