import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import is_BreedingValueMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import check_is_BreedingValueMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import cond_check_is_BreedingValueMatrix

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def mat():
    yield BreedingValueMatrix()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(BreedingValueMatrix)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(BreedingValueMatrix, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_location_is_abstract():
    generic_assert_abstract_property(BreedingValueMatrix, "location")

def test_scale_is_abstract():
    generic_assert_abstract_property(BreedingValueMatrix, "scale")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_targmax_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "targmax")

def test_targmin_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "targmin")

def test_tmax_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "tmax")

def test_tmean_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "tmean")

def test_tmin_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "tmin")

def test_trange_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "trange")

def test_tstd_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "tstd")

def test_tvar_is_abstract():
    generic_assert_abstract_method(BreedingValueMatrix, "tvar")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_BreedingValueMatrix_is_concrete():
    generic_assert_concrete_function(is_BreedingValueMatrix)

def test_check_is_BreedingValueMatrix_is_concrete():
    generic_assert_concrete_function(check_is_BreedingValueMatrix)

def test_cond_check_is_BreedingValueMatrix_is_concrete():
    generic_assert_concrete_function(cond_check_is_BreedingValueMatrix)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_BreedingValueMatrix(mat):
    assert is_BreedingValueMatrix(mat)

def test_check_is_BreedingValueMatrix(mat):
    with not_raises(TypeError):
        check_is_BreedingValueMatrix(mat, "mat")
    with pytest.raises(TypeError):
        check_is_BreedingValueMatrix(None, "mat")
