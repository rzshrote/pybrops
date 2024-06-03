import numpy
import pytest

from pybrops.core.mat.DenseScaledMatrix import DenseScaledMatrix, check_is_DenseScaledMatrix
from pybrops.test.assert_python import assert_class_isconcrete, assert_function_isconcrete, assert_method_isconcrete, assert_module_documentation, assert_module_public_api, assert_property_isconcrete, not_raises

################################################################################
################################ Test fixtures #################################
################################################################################

############################################################
###################### Scaled Matrix #######################
############################################################

@pytest.fixture
def length():
    yield 7

@pytest.fixture
def smat_mat(length):
    out = numpy.random.random((13,11,length))
    return out

@pytest.fixture
def smat_location(length):
    out = numpy.random.random((length,))
    yield out

@pytest.fixture
def smat_scale(length):
    out = numpy.random.random((length,))
    yield out

@pytest.fixture
def smat(
        smat_mat,
        smat_location,
        smat_scale,
    ):
    out = DenseScaledMatrix(
        mat = smat_mat,
        location = smat_location,
        scale = smat_scale,
    )
    yield out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_DenseScaledMatrix_module_documentation():
    import pybrops.core.mat.DenseScaledMatrix
    assert_module_documentation(pybrops.core.mat.DenseScaledMatrix)

def test_DenseScaledMatrix_module_public_api():
    import pybrops.core.mat.DenseScaledMatrix
    assert_module_public_api(pybrops.core.mat.DenseScaledMatrix)

################################################################################
############################ Test class attributes #############################
################################################################################

def test_DenseScaledMatrix_is_concrete():
    assert_class_isconcrete(DenseScaledMatrix)

################################################################################
######################## Test abstract special methods #########################
################################################################################

################################################################################
######################## Test concrete special methods #########################
################################################################################

### __init__
def test___init___is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "__init__")

def test___init__(smat_mat, smat_location, smat_scale):
    obj = DenseScaledMatrix(smat_mat, smat_location, smat_scale)
    assert isinstance(obj, DenseScaledMatrix)

### __copy__
def test___copy___is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "__copy__")

def test___copy__(smat):
    original = smat
    the_copy = smat.__copy__()
    assert id(original.mat) != id(the_copy.mat)
    assert id(original.location) != id(the_copy.location)
    assert id(original.scale) != id(the_copy.scale)
    assert numpy.all(original.mat == the_copy.mat)
    assert numpy.all(original.location == the_copy.location)
    assert numpy.all(original.scale == the_copy.scale)

### __deepcopy__
def test___deepcopy___is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "__deepcopy__")

def test___deepcopy__(smat):
    original = smat
    the_copy = smat.__deepcopy__()
    assert id(original.mat) != id(the_copy.mat)
    assert id(original.location) != id(the_copy.location)
    assert id(original.scale) != id(the_copy.scale)
    assert numpy.all(original.mat == the_copy.mat)
    assert numpy.all(original.location == the_copy.location)
    assert numpy.all(original.scale == the_copy.scale)

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
########################### Test concrete properties ###########################
################################################################################

### location
def test_location_is_concrete():
    assert_property_isconcrete(DenseScaledMatrix, "location")

### scale
def test_scale_is_concrete():
    assert_property_isconcrete(DenseScaledMatrix, "scale")

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
############################# Test concrete methods ############################
################################################################################

###################### Matrix copying ######################

### copy

def test_copy_is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "copy")

def test_copy(smat):
    original = smat
    the_copy = smat.copy()
    assert id(original.mat) != id(the_copy.mat)
    assert id(original.location) != id(the_copy.location)
    assert id(original.scale) != id(the_copy.scale)
    assert numpy.all(original.mat == the_copy.mat)
    assert numpy.all(original.location == the_copy.location)
    assert numpy.all(original.scale == the_copy.scale)

### deepcopy

def test_deepcopy_is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "deepcopy")

def test_deepcopy(smat):
    original = smat
    the_copy = smat.deepcopy()
    assert id(original.mat) != id(the_copy.mat)
    assert id(original.location) != id(the_copy.location)
    assert id(original.scale) != id(the_copy.scale)
    assert numpy.all(original.mat == the_copy.mat)
    assert numpy.all(original.location == the_copy.location)
    assert numpy.all(original.scale == the_copy.scale)

##################### Scaling methods ######################

### transform
def test_transform_is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "transform")

def test_transform(smat, smat_mat, smat_location, smat_scale):
    with not_raises(Exception):
        smat.transform(smat_mat)
    observed = smat.transform(smat_mat, copy = True)
    expected = (smat_mat - smat_location) / smat_scale
    assert numpy.allclose(observed, expected)

def test_transform_TypeError(smat, smat_mat):
    with pytest.raises(TypeError):
        smat.transform(object(), False)
    with pytest.raises(TypeError):
        smat.transform(smat_mat, object())

### untransform
def test_untransform_is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "untransform")

def test_untransform(smat, smat_mat, smat_location, smat_scale):
    with not_raises(Exception):
        smat.untransform(smat_mat)
    observed = smat.untransform(smat_mat, copy = True)
    expected = (smat_mat * smat_scale) + smat_location
    assert numpy.allclose(observed, expected)

def test_untransform_TypeError(smat, smat_mat):
    with pytest.raises(TypeError):
        smat.untransform(object(), False)
    with pytest.raises(TypeError):
        smat.untransform(smat_mat, object())

### rescale
def test_rescale_is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "rescale")

def test_rescale_not_inplace(smat, smat_mat):
    with not_raises(Exception):
        observed = smat.rescale(False)
    axes = tuple(range(observed.ndim-1))
    assert numpy.allclose(smat.mat, smat_mat)
    assert numpy.allclose(numpy.mean(observed, axis = axes), 0.0)
    assert numpy.allclose(numpy.std(observed, axis = axes),  1.0)

def test_rescale_inplace(smat, smat_mat):
    with not_raises(Exception):
        observed = smat.rescale(True)
    axes = tuple(range(observed.ndim-1))
    assert numpy.allclose(smat.mat, observed)
    assert numpy.allclose(numpy.mean(observed, axis = axes), 0.0)
    assert numpy.allclose(numpy.std(observed, axis = axes),  1.0)

def test_rescale_TypeError(smat):
    with pytest.raises(TypeError):
        smat.rescale(object())

### unscale
def test_unscale_is_concrete():
    assert_method_isconcrete(DenseScaledMatrix, "unscale")

def test_unscale_not_inplace(smat, smat_mat):
    with not_raises(Exception):
        observed = smat.unscale(False)
    axes = tuple(range(observed.ndim-1))
    assert numpy.allclose(smat.mat, smat_mat)
    assert numpy.all(observed != smat_mat)
    assert numpy.all(numpy.mean(observed, axis = axes) != 0.0)
    assert numpy.all(numpy.std(observed, axis = axes)  != 1.0)

def test_unscale_inplace(smat):
    smat_cpy = smat.mat.copy()
    with not_raises(Exception):
        observed = smat.unscale(True)
    axes = tuple(range(observed.ndim-1))
    assert numpy.all(smat.mat != smat_cpy)
    assert numpy.all(observed != smat_cpy)
    assert numpy.all(observed == smat.mat)
    assert numpy.all(numpy.mean(observed, axis = axes) != 0.0)
    assert numpy.all(numpy.std(observed, axis = axes)  != 1.0)

def test_unscale_TypeError(smat):
    with pytest.raises(TypeError):
        smat.unscale(object())

################################################################################
########################## Test abstract classmethods ##########################
################################################################################

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

### check_is_ScaledMatrix
def test_check_is_DenseScaledMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseScaledMatrix)

def test_check_is_DenseScaledMatrix(smat):
    with not_raises(TypeError):
        check_is_DenseScaledMatrix(smat, "smat")
    with pytest.raises(TypeError):
        check_is_DenseScaledMatrix(None, "smat")
