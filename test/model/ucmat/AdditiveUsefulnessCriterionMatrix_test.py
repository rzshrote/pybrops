import pytest

from pybrops.model.ucmat.AdditiveUsefulnessCriterionMatrix import AdditiveUsefulnessCriterionMatrix
from pybrops.model.ucmat.AdditiveUsefulnessCriterionMatrix import check_is_AdditiveUsefulnessCriterionMatrix
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
def ucmat():
    out = DummyAdditiveUsefulnessCriterionMatrix()
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_ProgenyMeanGenomicEstimatedBreedingValueMatrix_module_documentation():
    import pybrops.model.ucmat.AdditiveUsefulnessCriterionMatrix
    assert_module_documentation(pybrops.model.ucmat.AdditiveUsefulnessCriterionMatrix)

def test_ProgenyMeanGenomicEstimatedBreedingValueMatrix_module_public_api():
    import pybrops.model.ucmat.AdditiveUsefulnessCriterionMatrix
    assert_module_public_api(pybrops.model.ucmat.AdditiveUsefulnessCriterionMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_ProgenyMeanGenomicEstimatedBreedingValueMatrix_is_abstract():
    assert_class_isabstract(AdditiveUsefulnessCriterionMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

### epgc
def test_epgc_is_abstract():
    assert_property_isabstract(AdditiveUsefulnessCriterionMatrix, "epgc")

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
    assert_classmethod_isabstract(AdditiveUsefulnessCriterionMatrix, "from_gmod")

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
    assert_function_isconcrete(check_is_AdditiveUsefulnessCriterionMatrix)

def test_check_is_ProgenyMeanGenomicEstimatedBreedingValueMatrix(ucmat):
    with not_raises(TypeError):
        check_is_AdditiveUsefulnessCriterionMatrix(ucmat, "ucmat")
    with pytest.raises(TypeError):
        check_is_AdditiveUsefulnessCriterionMatrix(None, "ucmat")