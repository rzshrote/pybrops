import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

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
    assert_docstring(SquareTaxaMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SquareTaxaMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nsquare_taxa_is_abstract():
    assert_abstract_property(SquareTaxaMatrix, "nsquare_taxa")

def test_square_taxa_axes_is_abstract():
    assert_abstract_property(SquareTaxaMatrix, "square_taxa_axes")

def test_square_taxa_axes_len_is_abstract():
    assert_abstract_property(SquareTaxaMatrix, "square_taxa_axes_len")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_is_square_taxa_is_abstract():
    assert_abstract_method(SquareTaxaMatrix, "is_square_taxa")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SquareTaxaMatrix_is_concrete():
    assert_concrete_function(check_is_SquareTaxaMatrix)

def test_check_is_SquareTaxaMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareTaxaMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareTaxaMatrix(None, "mat")