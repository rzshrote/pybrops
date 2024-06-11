import pytest

from pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix import ParentalMeanEstimatedBreedingValueMatrix
from pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix import check_is_ParentalMeanEstimatedBreedingValueMatrix
from pybrops.test.assert_python import assert_class_isabstract
from pybrops.test.assert_python import assert_classmethod_isabstract
from pybrops.test.assert_python import assert_function_isconcrete
from pybrops.test.assert_python import assert_module_documentation
from pybrops.test.assert_python import assert_module_public_api
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import not_raises

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
    import pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix
    assert_module_documentation(pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix)

def test_ProgenyMeanEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix
    assert_module_public_api(pybrops.model.pmebvmat.ParentalMeanEstimatedBreedingValueMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_ProgenyMeanEstimatedBreedingValueMatrix_is_abstract():
    assert_class_isabstract(ParentalMeanEstimatedBreedingValueMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

### epgc
def test_epgc_is_abstract():
    assert_property_isabstract(ParentalMeanEstimatedBreedingValueMatrix, "epgc")

################################################################################
######################## Test concrete special methods #########################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

### from_bvmat
def test_from_bvmat_is_abstract():
    assert_classmethod_isabstract(ParentalMeanEstimatedBreedingValueMatrix, "from_bvmat")

################################################################################
########################## Test concrete classmethods ##########################
################################################################################

################################################################################
######################### Test abstract staticmethods ##########################
################################################################################

################################################################################
######################### Test concrete staticmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_ProgenyMeanEstimatedBreedingValueMatrix
def test_check_is_ProgenyMeanEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_ParentalMeanEstimatedBreedingValueMatrix)

def test_check_is_ProgenyMeanEstimatedBreedingValueMatrix(pmebvmat):
    with not_raises(TypeError):
        check_is_ParentalMeanEstimatedBreedingValueMatrix(pmebvmat, "pmebvmat")
    with pytest.raises(TypeError):
        check_is_ParentalMeanEstimatedBreedingValueMatrix(None, "pmebvmat")
