from pybrops.core.mat.ScaledMatrix import check_is_ScaledMatrix
from pybrops.test.assert_python import assert_class_isabstract, assert_function_isconcrete, assert_method_isabstract, assert_module_documentation, assert_module_public_api, assert_property_isabstract, not_raises
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
###################### Scaled Matrix #######################
############################################################

@pytest.fixture
def smat():
    out = DummyScaledMatrix()
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_ScaledMatrix_module_documentation():
    import pybrops.core.mat.ScaledMatrix
    assert_module_documentation(pybrops.core.mat.ScaledMatrix)

def test_ScaledMatrix_module_public_api():
    import pybrops.core.mat.ScaledMatrix
    assert_module_public_api(pybrops.core.mat.ScaledMatrix)

################################################################################
########################### Test class documentation ###########################
################################################################################

def test_ScaledMatrix_is_abstract():
    assert_class_isabstract(ScaledMatrix)

################################################################################
############################ Test Class Properties #############################
################################################################################

### location
def test_location_is_abstract():
    assert_property_isabstract(ScaledMatrix, "location")

### scale
def test_scale_is_abstract():
    assert_property_isabstract(ScaledMatrix, "scale")

################################################################################
########################## Test Class Special Methods ##########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

### transform
def test_transform_is_abstract():
    assert_method_isabstract(ScaledMatrix, "transform")

### untransform
def test_untransform_is_abstract():
    assert_method_isabstract(ScaledMatrix, "untransform")

### rescale
def test_rescale_is_abstract():
    assert_method_isabstract(ScaledMatrix, "rescale")

### unscale
def test_unscale_is_abstract():
    assert_method_isabstract(ScaledMatrix, "unscale")

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################

### check_is_ScaledMatrix
def test_check_is_ScaledMatrix_is_concrete():
    assert_function_isconcrete(check_is_ScaledMatrix)

def test_check_is_ScaledMatrix(smat):
    with not_raises(TypeError):
        check_is_ScaledMatrix(smat, "smat")
    with pytest.raises(TypeError):
        check_is_ScaledMatrix(None, "smat")
