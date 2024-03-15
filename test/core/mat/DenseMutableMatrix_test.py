import pytest
import numpy

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybrops.core.mat.DenseMutableMatrix import check_is_DenseMutableMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat_float64():
    a = numpy.float64([
        [3.3, 9.2, 5.6],
        [8.7, 3.7, 4.1],
        [9.0, 4.7, 3.8]
    ])
    yield a

@pytest.fixture
def mat(mat_float64):
    yield DenseMutableMatrix(mat_float64)

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(DenseMutableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test___init___is_concrete():
    assert_method_isconcrete(DenseMutableMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_DenseMutableMatrix_is_concrete():
    assert_function_isconcrete(check_is_DenseMutableMatrix)

def test_check_is_DenseMutableMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseMutableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseMutableMatrix(None, "mat")
