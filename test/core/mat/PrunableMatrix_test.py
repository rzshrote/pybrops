import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.PrunableMatrix import PrunableMatrix
from pybrops.core.mat.PrunableMatrix import check_is_PrunableMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyPrunableMatrix()

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
def test_prune_is_abstract():
    assert_abstract_method(PrunableMatrix, "prune")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PrunableMatrix_is_concrete():
    assert_concrete_function(check_is_PrunableMatrix)

def test_check_is_PrunableMatrix(mat):
    with not_raises(TypeError):
        check_is_PrunableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PrunableMatrix(None, "mat")
