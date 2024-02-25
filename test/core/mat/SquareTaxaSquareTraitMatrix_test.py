import pytest
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete
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
    assert_class_documentation(SquareTaxaSquareTraitMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

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
    assert_function_isconcrete(check_is_SquareTaxaSquareTraitMatrix)

def test_check_is_SquareTaxaSquareTraitMatrix(mat):
    with not_raises(TypeError):
        check_is_SquareTaxaSquareTraitMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_SquareTaxaSquareTraitMatrix(None, "mat")
