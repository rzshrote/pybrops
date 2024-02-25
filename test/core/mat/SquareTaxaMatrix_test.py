import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.SquareTaxaMatrix import check_is_SquareTaxaMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummySquareTaxaMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(SquareTaxaMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nsquare_taxa_is_abstract():
    assert_property_isabstract(SquareTaxaMatrix, "nsquare_taxa")

def test_square_taxa_axes_is_abstract():
    assert_property_isabstract(SquareTaxaMatrix, "square_taxa_axes")

def test_square_taxa_axes_len_is_abstract():
    assert_property_isabstract(SquareTaxaMatrix, "square_taxa_axes_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_is_square_taxa_is_abstract():
    assert_method_isabstract(SquareTaxaMatrix, "is_square_taxa")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SquareTaxaMatrix_is_concrete():
    assert_function_isconcrete(check_is_SquareTaxaMatrix)

def test_check_is_SquareTaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareTaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareTaxaMatrix(None, "mat")
