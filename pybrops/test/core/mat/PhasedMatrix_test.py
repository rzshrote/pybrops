import inspect
import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.core.mat.PhasedMatrix import PhasedMatrix
from pybrops.core.mat.PhasedMatrix import check_is_PhasedMatrix
from pybrops.test.core.mat.common_fixtures import *

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
    assert_docstring(PhasedMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(PhasedMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nphase_is_abstract():
    assert_abstract_property(PhasedMatrix, "nphase")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "adjoin_phase")

def test_delete_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "delete_phase")

def test_insert_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "insert_phase")

def test_select_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "select_phase")

def test_concat_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "concat_phase")

def test_append_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "append_phase")

def test_remove_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "remove_phase")

def test_incorp_phase_is_abstract():
    assert_abstract_method(PhasedMatrix, "incorp_phase")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_PhasedMatrix_is_concrete():
    assert_concrete_function(check_is_PhasedMatrix)

def test_check_is_PhasedMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedMatrix(None, "mat")
