import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import is_GeneticVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import check_is_GeneticVarianceMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield GeneticVarianceMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(GeneticVarianceMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GeneticVarianceMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_append_is_abstract(mat):
    assert_abstract_method(mat, "from_gmod")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GeneticVarianceMatrix_is_concrete():
    assert_concrete_function(is_GeneticVarianceMatrix)

def test_is_GeneticVarianceMatrix(mat):
    assert is_GeneticVarianceMatrix(mat)

def test_check_is_GeneticVarianceMatrix_is_concrete():
    assert_concrete_function(check_is_GeneticVarianceMatrix)

def test_check_is_GeneticVarianceMatrix(mat):
    with not_raises(TypeError):
        check_is_GeneticVarianceMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GeneticVarianceMatrix(None, "mat")
