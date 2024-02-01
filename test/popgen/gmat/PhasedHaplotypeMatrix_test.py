import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.popgen.gmat.PhasedHaplotypeMatrix import PhasedHaplotypeMatrix
from pybrops.popgen.gmat.PhasedHaplotypeMatrix import check_is_PhasedHaplotypeMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyPhasedHaplotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(PhasedHaplotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(PhasedHaplotypeMatrix, "__init__")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_PhasedHaplotypeMatrix_is_concrete():
    assert_concrete_function(check_is_PhasedHaplotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PhasedHaplotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedHaplotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedHaplotypeMatrix(None, "mat")
