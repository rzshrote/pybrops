import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.popgen.gmat import PhasedHaplotypeMatrix
from pybropt.popgen.gmat import is_PhasedHaplotypeMatrix
from pybropt.popgen.gmat import check_is_PhasedHaplotypeMatrix
from pybropt.popgen.gmat import cond_check_is_PhasedHaplotypeMatrix


################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield PhasedHaplotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(PhasedHaplotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(PhasedHaplotypeMatrix, "__init__")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_PhasedHaplotypeMatrix_is_concrete():
    generic_assert_concrete_function(is_PhasedHaplotypeMatrix)

def test_check_is_PhasedHaplotypeMatrix_is_concrete():
    generic_assert_concrete_function(check_is_PhasedHaplotypeMatrix)

def test_cond_check_is_PhasedHaplotypeMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_PhasedHaplotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhasedHaplotypeMatrix(mat):
    assert is_PhasedHaplotypeMatrix(mat)

def test_check_is_PhasedHaplotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedHaplotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedHaplotypeMatrix(None, "mat")
