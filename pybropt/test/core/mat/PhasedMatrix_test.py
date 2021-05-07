import inspect
import pytest

from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.core.mat import PhasedMatrix
from pybropt.core.mat import is_PhasedMatrix
from pybropt.core.mat import check_is_PhasedMatrix
from pybropt.core.mat import cond_check_is_PhasedMatrix

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
    generic_assert_docstring(PhasedMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(PhasedMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_nphase_is_abstract():
    generic_assert_abstract_property(PhasedMatrix, "nphase")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_adjoin_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "adjoin_phase")

def test_delete_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "delete_phase")

def test_insert_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "insert_phase")

def test_select_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "select_phase")

def test_concat_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "concat_phase")

def test_append_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "append_phase")

def test_remove_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "remove_phase")

def test_incorp_phase_is_abstract(mat):
    generic_assert_abstract_method(mat, "incorp_phase")

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhasedMatrix_is_concrete():
    generic_assert_concrete_function(is_PhasedMatrix)

def test_is_PhasedMatrix(mat):
    assert is_PhasedMatrix(mat)

def test_check_is_PhasedMatrix_is_concrete():
    generic_assert_concrete_function(check_is_PhasedMatrix)

def test_check_is_PhasedMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedMatrix(None, "mat")

def test_cond_check_is_PhasedMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_PhasedMatrix)
