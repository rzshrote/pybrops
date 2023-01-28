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

from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import is_PhasedGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import check_is_PhasedGenotypeMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield PhasedGenotypeMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(PhasedGenotypeMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(PhasedGenotypeMatrix, "__init__")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_PhasedGenotypeMatrix_is_concrete():
    generic_assert_concrete_function(is_PhasedGenotypeMatrix)

def test_check_is_PhasedGenotypeMatrix_is_concrete():
    generic_assert_concrete_function(check_is_PhasedGenotypeMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_PhasedGenotypeMatrix(mat):
    assert is_PhasedGenotypeMatrix(mat)

def test_check_is_PhasedGenotypeMatrix(mat):
    with not_raises(TypeError):
        check_is_PhasedGenotypeMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_PhasedGenotypeMatrix(None, "mat")
