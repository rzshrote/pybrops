import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function
from .common_fixtures import *
from pybrops.core.mat.SquareTaxaSquareTraitMatrix import SquareTaxaSquareTraitMatrix
from pybrops.core.mat.SquareTaxaSquareTraitMatrix import check_is_SquareTaxaSquareTraitMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummySquareTaxaSquareTraitMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(SquareTaxaSquareTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(SquareTaxaSquareTraitMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_SquareTaxaSquareTraitMatrix_is_concrete():
    assert_concrete_function(check_is_SquareTaxaSquareTraitMatrix)

def test_check_is_SquareTaxaSquareTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareTaxaSquareTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareTaxaSquareTraitMatrix(None, "mat")
