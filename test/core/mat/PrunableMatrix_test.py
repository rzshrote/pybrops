import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

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
    assert_class_documentation(PrunableMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_prune_is_abstract():
    assert_method_isabstract(PrunableMatrix, "prune")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PrunableMatrix_is_concrete():
    assert_function_isconcrete(check_is_PrunableMatrix)

def test_check_is_PrunableMatrix(mat):
    with not_raises(TypeError):
        check_is_PrunableMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PrunableMatrix(None, "mat")
