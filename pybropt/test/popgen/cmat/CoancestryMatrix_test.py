import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.popgen.cmat import CoancestryMatrix
from pybropt.popgen.cmat import is_CoancestryMatrix
from pybropt.popgen.cmat import check_is_CoancestryMatrix
from pybropt.popgen.cmat import cond_check_is_CoancestryMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def cmat():
    yield CoancestryMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(CoancestryMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(CoancestryMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_coancestry_is_abstract():
    generic_assert_abstract_method(CoancestryMatrix, "coancestry")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_CoancestryMatrix_is_concrete():
    generic_assert_concrete_function(is_CoancestryMatrix)

def test_check_is_CoancestryMatrix_is_concrete():
    generic_assert_concrete_function(check_is_CoancestryMatrix)

def test_cond_check_is_CoancestryMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_CoancestryMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_CoancestryMatrix(cmat):
    assert is_CoancestryMatrix(cmat)

def test_check_is_CoancestryMatrix(cmat):
    with not_raises(TypeError):
        check_is_CoancestryMatrix(cmat, "cmat")
    with pytest.raises(TypeError):
        check_is_CoancestryMatrix(None, "cmat")
