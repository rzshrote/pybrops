import pytest

from pybrops.model.pmgebvmat.ParentalMeanGenomicEstimatedBreedingValueMatrix import ParentalMeanGenomicEstimatedBreedingValueMatrix
from pybrops.model.pmgebvmat.ParentalMeanGenomicEstimatedBreedingValueMatrix import check_is_ParentalMeanGenomicEstimatedBreedingValueMatrix
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
def pmgebvmat():
    out = DummyProgenyMeanGenomicEstimatedBreedingValueMatrix()
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_ProgenyMeanGenomicEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.pmgebvmat.ParentalMeanGenomicEstimatedBreedingValueMatrix
    assert_module_documentation(pybrops.model.pmgebvmat.ParentalMeanGenomicEstimatedBreedingValueMatrix)

def test_ProgenyMeanGenomicEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.pmgebvmat.ParentalMeanGenomicEstimatedBreedingValueMatrix
    assert_module_public_api(pybrops.model.pmgebvmat.ParentalMeanGenomicEstimatedBreedingValueMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_ProgenyMeanGenomicEstimatedBreedingValueMatrix_is_abstract():
    assert_class_isabstract(ParentalMeanGenomicEstimatedBreedingValueMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

### epgc
def test_epgc_is_abstract():
    assert_property_isabstract(ParentalMeanGenomicEstimatedBreedingValueMatrix, "epgc")

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

### from_gmod
def test_from_gmod_is_abstract():
    assert_classmethod_isabstract(ParentalMeanGenomicEstimatedBreedingValueMatrix, "from_gmod")

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

### check_is_ProgenyMeanGenomicEstimatedBreedingValueMatrix
def test_check_is_ProgenyMeanGenomicEstimatedBreedingValueMatrix_is_concrete():
    assert_function_isconcrete(check_is_ParentalMeanGenomicEstimatedBreedingValueMatrix)

def test_check_is_ProgenyMeanGenomicEstimatedBreedingValueMatrix(pmgebvmat):
    with not_raises(TypeError):
        check_is_ParentalMeanGenomicEstimatedBreedingValueMatrix(pmgebvmat, "pmgebvmat")
    with pytest.raises(TypeError):
        check_is_ParentalMeanGenomicEstimatedBreedingValueMatrix(None, "pmgebvmat")
