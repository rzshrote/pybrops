import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.breed.prot.bv.BreedingValueProtocol import is_BreedingValueProtocol
from pybrops.breed.prot.bv.BreedingValueProtocol import check_is_BreedingValueProtocol

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def bvprot():
    yield BreedingValueProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(BreedingValueProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(BreedingValueProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_estimate_is_abstract():
    generic_assert_abstract_method(BreedingValueProtocol, "estimate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_BreedingValueProtocol_is_concrete():
    generic_assert_concrete_function(is_BreedingValueProtocol)

def test_check_is_BreedingValueProtocol_is_concrete():
    generic_assert_concrete_function(check_is_BreedingValueProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_BreedingValueProtocol(bvprot):
    assert is_BreedingValueProtocol(bvprot)

def test_check_is_BreedingValueProtocol(bvprot):
    with not_raises(TypeError):
        check_is_BreedingValueProtocol(bvprot, "bvprot")
    with pytest.raises(TypeError):
        check_is_BreedingValueProtocol(None, "bvprot")
