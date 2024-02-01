import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

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
    assert_docstring(AdditiveGenicVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(AdditiveGenicVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract():
    assert_abstract_method(AdditiveGenicVarianceMatrix, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_AdditiveGenicVarianceMatrix_is_concrete():
    assert_concrete_function(check_is_AdditiveGenicVarianceMatrix)

def test_check_is_AdditiveGenicVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_AdditiveGenicVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_AdditiveGenicVarianceMatrix(None, "mat")
