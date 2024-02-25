import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import check_is_AdditiveGenicVarianceMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyAdditiveGenicVarianceMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(AdditiveGenicVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(AdditiveGenicVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract():
    assert_method_isabstract(AdditiveGenicVarianceMatrix, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveGenicVarianceMatrix_is_concrete():
    assert_function_isconcrete(check_is_AdditiveGenicVarianceMatrix)

def test_check_is_AdditiveGenicVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_AdditiveGenicVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_AdditiveGenicVarianceMatrix(None, "mat")
