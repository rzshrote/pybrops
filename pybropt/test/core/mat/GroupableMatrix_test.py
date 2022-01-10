import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat.GroupableMatrix import GroupableMatrix
from pybropt.core.mat.GroupableMatrix import is_GroupableMatrix
from pybropt.core.mat.GroupableMatrix import check_is_GroupableMatrix
from pybropt.core.mat.GroupableMatrix import cond_check_is_GroupableMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield GroupableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GroupableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GroupableMatrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_group_is_abstract(mat):
    generic_assert_abstract_method(mat, "group")

def test_is_grouped_is_abstract(mat):
    generic_assert_abstract_method(mat, "is_grouped")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GroupableMatrix_is_concrete():
    generic_assert_concrete_function(is_GroupableMatrix)

def test_is_GroupableMatrix(mat):
    assert is_GroupableMatrix(mat)

def test_check_is_GroupableMatrix_is_concrete():
    generic_assert_concrete_function(check_is_GroupableMatrix)

def test_check_is_GroupableMatrix(mat):
    with not_raises(TypeError):
        check_is_GroupableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GroupableMatrix(None, "mat")

def test_cond_check_is_GroupableMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_GroupableMatrix)
