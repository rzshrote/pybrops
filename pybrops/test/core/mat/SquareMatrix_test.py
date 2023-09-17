import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.SquareMatrix import check_is_SquareMatrix
from pybrops.test.core.mat.common_fixtures import *

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
    assert_docstring(SquareMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SquareMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_mat_is_abstract():
    assert_abstract_property(SquareMatrix, "nsquare")

def test_mat_is_abstract():
    assert_abstract_property(SquareMatrix, "square_axes")

def test_mat_is_abstract():
    assert_abstract_property(SquareMatrix, "square_axes_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract():
    assert_abstract_method(SquareMatrix, "is_square")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SquareMatrix_is_concrete():
    assert_concrete_function(check_is_SquareMatrix)

def test_check_is_SquareMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareMatrix(None, "mat")
