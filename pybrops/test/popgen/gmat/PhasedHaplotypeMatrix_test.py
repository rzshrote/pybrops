import inspect
import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.gmat.PhasedHaplotypeMatrix import PhasedHaplotypeMatrix
from pybrops.popgen.gmat.PhasedHaplotypeMatrix import is_PhasedHaplotypeMatrix
from pybrops.popgen.gmat.PhasedHaplotypeMatrix import check_is_PhasedHaplotypeMatrix

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
