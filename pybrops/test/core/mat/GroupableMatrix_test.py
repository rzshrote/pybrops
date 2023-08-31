import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.GroupableMatrix import GroupableMatrix
from pybrops.core.mat.GroupableMatrix import check_is_GroupableMatrix
from pybrops.test.core.mat.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyGroupableMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(GroupableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GroupableMatrix, "__init__")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_group_is_abstract():
    assert_abstract_method(GroupableMatrix, "group")

def test_is_grouped_is_abstract():
    assert_abstract_method(GroupableMatrix, "is_grouped")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GroupableMatrix_is_concrete():
    assert_concrete_function(check_is_GroupableMatrix)

def test_check_is_GroupableMatrix(mat):
    with not_raises(TypeError):
        check_is_GroupableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GroupableMatrix(None, "mat")
