import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.GroupableMatrix import GroupableMatrix
from pybrops.core.mat.GroupableMatrix import check_is_GroupableMatrix
from .common_fixtures import *

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
    assert_class_documentation(GroupableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_group_is_abstract():
    assert_method_isabstract(GroupableMatrix, "group")

def test_ungroup_is_abstract():
    assert_method_isabstract(GroupableMatrix, "ungroup")

def test_is_grouped_is_abstract():
    assert_method_isabstract(GroupableMatrix, "is_grouped")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GroupableMatrix_is_concrete():
    assert_function_isconcrete(check_is_GroupableMatrix)

def test_check_is_GroupableMatrix(mat):
    with not_raises(TypeError):
        check_is_GroupableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_GroupableMatrix(None, "mat")
