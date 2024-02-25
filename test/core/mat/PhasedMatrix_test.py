import inspect
import pytest

from pybrops.test.assert_python import assert_classmethod_isabstract, not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_function_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.core.mat.PhasedMatrix import PhasedMatrix
from pybrops.core.mat.PhasedMatrix import check_is_PhasedMatrix
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield DummyPhasedMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(PhasedMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nphase_is_abstract():
    assert_property_isabstract(PhasedMatrix, "nphase")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "adjoin_phase")

def test_delete_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "delete_phase")

def test_insert_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "insert_phase")

def test_select_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "select_phase")

def test_concat_phase_is_abstract():
    assert_classmethod_isabstract(PhasedMatrix, "concat_phase")

def test_append_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "append_phase")

def test_remove_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "remove_phase")

def test_incorp_phase_is_abstract():
    assert_method_isabstract(PhasedMatrix, "incorp_phase")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PhasedMatrix_is_concrete():
    assert_function_isconcrete(check_is_PhasedMatrix)

def test_check_is_PhasedMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedMatrix(None, "mat")
