import numpy
import pytest

from pybrops.core.error.error_value_numpy import check_ndarray_is_square_along_axes
from pybrops.test.assert_python import assert_function_isconcrete, assert_module_documentation, assert_module_public_api

################################################################################
################################ Test fixtures #################################
################################################################################

@pytest.fixture
def ndarray_nonsquare():
    out = numpy.random.random((2,3,5,7))
    return out

@pytest.fixture
def ndarray_square():
    out = numpy.random.random((4,4,4,4))
    return out

@pytest.fixture
def ndarray_square2():
    out = numpy.random.random((2,2,3,3))
    return out

################################################################################
########################## Test module documentation ###########################
################################################################################

def test_error_value_numpy_module_documentation():
    import pybrops.core.error.error_value_numpy
    assert_module_documentation(pybrops.core.error.error_value_numpy)

def test_error_value_numpy_module_public_api():
    import pybrops.core.error.error_value_numpy
    assert_module_public_api(pybrops.core.error.error_value_numpy)

################################################################################
############################ Test module functions #############################
################################################################################

### check_ndarray_is_square_along_axes

def test_check_ndarray_is_square_along_axes_is_concrete():
    assert_function_isconcrete(check_ndarray_is_square_along_axes)

def test_check_ndarray_is_square_along_axes_ValueError(
        ndarray_nonsquare, 
        ndarray_square,
        ndarray_square2,
    ):
    with pytest.raises(ValueError):
        check_ndarray_is_square_along_axes(ndarray_square, "square", (0,))
    with pytest.raises(ValueError):
        check_ndarray_is_square_along_axes(ndarray_nonsquare, "nonsquare")
    with pytest.raises(ValueError):
        check_ndarray_is_square_along_axes(ndarray_square2, "square2", (1,2))

def test_check_ndarray_is_square_along_axes(
        ndarray_square,
        ndarray_square2,
    ):
    assert check_ndarray_is_square_along_axes(ndarray_square, "square") is None
    assert check_ndarray_is_square_along_axes(ndarray_square, "square", (0,1)) is None
    assert check_ndarray_is_square_along_axes(ndarray_square, "square", (-1,-2)) is None
    assert check_ndarray_is_square_along_axes(ndarray_square2, "square2", (0,1)) is None
    assert check_ndarray_is_square_along_axes(ndarray_square2, "square2", (2,3)) is None
