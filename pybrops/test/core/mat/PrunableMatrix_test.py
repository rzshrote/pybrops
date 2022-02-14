import inspect
import pytest

from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.core.mat.PrunableMatrix import PrunableMatrix
from pybrops.core.mat.PrunableMatrix import is_PrunableMatrix
from pybrops.core.mat.PrunableMatrix import check_is_PrunableMatrix
from pybrops.core.mat.PrunableMatrix import cond_check_is_PrunableMatrix

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
    generic_assert_docstring(PrunableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(PrunableMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_prune_is_abstract(mat):
    generic_assert_abstract_method(mat, "prune")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PrunableMatrix_is_concrete():
    generic_assert_concrete_function(is_PrunableMatrix)

def test_is_PrunableMatrix(mat):
    assert is_PrunableMatrix(mat)

def test_check_is_PrunableMatrix_is_concrete():
    generic_assert_concrete_function(check_is_PrunableMatrix)

def test_check_is_PrunableMatrix(mat):
    with not_raises(TypeError):
        check_is_PrunableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PrunableMatrix(None, "mat")

def test_cond_check_is_PrunableMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_PrunableMatrix)

# def test_cond_check_is_PrunableMatrix():
#     with not_raises(TypeError):
#         cond_check_is_PrunableMatrix(mat, "mat")
#     # with not_raises(TypeError):
#     #     cond_check_is_PrunableMatrix(None, "mat")
#     with pytest.raises(TypeError):
#         cond_check_is_PrunableMatrix("string", "mat")
