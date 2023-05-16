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
from pybrops.core.mat.PhasedMatrix import is_PhasedMatrix
from pybrops.core.mat.PhasedMatrix import check_is_PhasedMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield PhasedMatrix()

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
def test_adjoin_phase_is_abstract(mat):
    assert_abstract_method(mat, "adjoin_phase")

def test_delete_phase_is_abstract(mat):
    assert_abstract_method(mat, "delete_phase")

def test_insert_phase_is_abstract(mat):
    assert_abstract_method(mat, "insert_phase")

def test_select_phase_is_abstract(mat):
    assert_abstract_method(mat, "select_phase")

def test_concat_phase_is_abstract(mat):
    assert_abstract_method(mat, "concat_phase")

def test_append_phase_is_abstract(mat):
    assert_abstract_method(mat, "append_phase")

def test_remove_phase_is_abstract(mat):
    assert_abstract_method(mat, "remove_phase")

def test_incorp_phase_is_abstract(mat):
    assert_abstract_method(mat, "incorp_phase")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhasedMatrix_is_concrete():
    assert_concrete_function(is_PhasedMatrix)

def test_is_PhasedMatrix(mat):
    assert is_PhasedMatrix(mat)

def test_check_is_PhasedMatrix_is_concrete():
    assert_concrete_function(check_is_PhasedMatrix)

def test_check_is_PhasedMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedMatrix(None, "mat")
