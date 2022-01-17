import pytest
import numpy

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybropt.core.mat.DenseMutableMatrix import is_DenseMutableMatrix
from pybropt.core.mat.DenseMutableMatrix import check_is_DenseMutableMatrix
from pybropt.core.mat.DenseMutableMatrix import cond_check_is_DenseMutableMatrix

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
    generic_assert_docstring(DenseMutableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(DenseMutableMatrix, "__init__")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_DenseMutableMatrix_is_concrete():
    generic_assert_concrete_function(is_DenseMutableMatrix)

def test_is_DenseMutableMatrix(mat):
    assert is_DenseMutableMatrix(mat)

def test_check_is_DenseMutableMatrix_is_concrete():
    generic_assert_concrete_function(check_is_DenseMutableMatrix)

def test_check_is_DenseMutableMatrix(mat):
    with not_raises(TypeError):
        check_is_DenseMutableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_DenseMutableMatrix(None, "mat")

def test_cond_check_is_DenseMutableMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_DenseMutableMatrix)
