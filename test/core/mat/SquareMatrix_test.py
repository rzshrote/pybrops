import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
from .common_fixtures import *
from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.SquareMatrix import check_is_SquareMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummySquareMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(SquareMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nsquare_is_abstract():
    assert_property_isabstract(SquareMatrix, "nsquare")

def test_square_axes_is_abstract():
    assert_property_isabstract(SquareMatrix, "square_axes")

def test_square_axes_len_is_abstract():
    assert_property_isabstract(SquareMatrix, "square_axes_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_is_square_is_abstract():
    assert_method_isabstract(SquareMatrix, "is_square")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SquareMatrix_is_concrete():
    assert_function_isconcrete(check_is_SquareMatrix)

def test_check_is_SquareMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareMatrix(None, "mat")
