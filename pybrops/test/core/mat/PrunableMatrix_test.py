import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.PrunableMatrix import PrunableMatrix
from pybrops.core.mat.PrunableMatrix import is_PrunableMatrix
from pybrops.core.mat.PrunableMatrix import check_is_PrunableMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield PrunableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(PrunableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(PrunableMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_prune_is_abstract(mat):
    assert_abstract_method(mat, "prune")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PrunableMatrix_is_concrete():
    assert_concrete_function(is_PrunableMatrix)

def test_is_PrunableMatrix(mat):
    assert is_PrunableMatrix(mat)

def test_check_is_PrunableMatrix_is_concrete():
    assert_concrete_function(check_is_PrunableMatrix)

def test_check_is_PrunableMatrix(mat):
    with not_raises(TypeError):
        check_is_PrunableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PrunableMatrix(None, "mat")
