import pytest
import numpy

from pybrops.test import not_raises
from pybrops.test import assert_docstring
from pybrops.test import assert_abstract_method
from pybrops.test import assert_abstract_function
from pybrops.test import assert_abstract_property
from pybrops.test import assert_concrete_method
from pybrops.test import assert_concrete_function

from pybrops.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybrops.core.mat.DenseMutableMatrix import is_DenseMutableMatrix
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
    assert_docstring(DenseMutableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(DenseMutableMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseMutableMatrix_is_concrete():
    assert_concrete_function(is_DenseMutableMatrix)

def test_is_DenseMutableMatrix(mat):
    assert is_DenseMutableMatrix(mat)

def test_check_is_DenseMutableMatrix_is_concrete():
    assert_concrete_function(check_is_DenseMutableMatrix)

def test_check_is_DenseMutableMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseMutableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseMutableMatrix(None, "mat")
