import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
from .common_fixtures import *
from pybrops.core.mat.SquareTraitMatrix import SquareTraitMatrix
from pybrops.core.mat.SquareTraitMatrix import check_is_SquareTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummySquareTraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(SquareTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nsquare_trait_is_abstract():
    assert_property_isabstract(SquareTraitMatrix, "nsquare_trait")

def test_square_trait_axes_is_abstract():
    assert_property_isabstract(SquareTraitMatrix, "square_trait_axes")

def test_square_trait_axes_len_is_abstract():
    assert_property_isabstract(SquareTraitMatrix, "square_trait_axes_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_is_square_trait_is_abstract():
    assert_method_isabstract(SquareTraitMatrix, "is_square_trait")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SquareTraitMatrix_is_concrete():
    assert_function_isconcrete(check_is_SquareTraitMatrix)

def test_check_is_SquareTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareTraitMatrix(None, "mat")
