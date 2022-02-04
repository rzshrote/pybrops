import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import is_AdditiveGeneticVarianceMatrix
from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import check_is_AdditiveGeneticVarianceMatrix
from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import cond_check_is_AdditiveGeneticVarianceMatrix

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
    generic_assert_docstring(AdditiveGeneticVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(AdditiveGeneticVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract(mat):
    generic_assert_abstract_method(mat, "from_algmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_AdditiveGeneticVarianceMatrix_is_concrete():
    generic_assert_concrete_function(is_AdditiveGeneticVarianceMatrix)

def test_is_AdditiveGeneticVarianceMatrix(mat):
    assert is_AdditiveGeneticVarianceMatrix(mat)

def test_check_is_AdditiveGeneticVarianceMatrix_is_concrete():
    generic_assert_concrete_function(check_is_AdditiveGeneticVarianceMatrix)

def test_check_is_AdditiveGeneticVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_AdditiveGeneticVarianceMatrix(None, "mat")

def test_cond_check_is_AdditiveGeneticVarianceMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_AdditiveGeneticVarianceMatrix)
