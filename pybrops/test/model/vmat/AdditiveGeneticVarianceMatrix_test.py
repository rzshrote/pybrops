import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import is_AdditiveGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import check_is_AdditiveGeneticVarianceMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield AdditiveGeneticVarianceMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(AdditiveGeneticVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(AdditiveGeneticVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_from_algmod_is_abstract(mat):
    assert_abstract_method(mat, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_AdditiveGeneticVarianceMatrix_is_concrete():
    assert_concrete_function(is_AdditiveGeneticVarianceMatrix)

def test_is_AdditiveGeneticVarianceMatrix(mat):
    assert is_AdditiveGeneticVarianceMatrix(mat)

def test_check_is_AdditiveGeneticVarianceMatrix_is_concrete():
    assert_concrete_function(check_is_AdditiveGeneticVarianceMatrix)

def test_check_is_AdditiveGeneticVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrix(None, "mat")
