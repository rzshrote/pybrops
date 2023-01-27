import pytest

from pybrops.test import generic_test_abstract_methods
from pybrops.test import not_raises
from pybrops.test import generic_assert_docstring
from pybrops.test import generic_assert_abstract_method
from pybrops.test import generic_assert_abstract_function
from pybrops.test import generic_assert_abstract_property
from pybrops.test import generic_assert_concrete_method
from pybrops.test import generic_assert_concrete_function

from pybrops.breed.arch.BreedingEdge import BreedingEdge
from pybrops.breed.arch.BreedingEdge import is_BreedingEdge
from pybrops.breed.arch.BreedingEdge import check_is_BreedingEdge

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield BreedingEdge()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(BreedingEdge)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(BreedingEdge, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(BreedingEdge, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_BreedingEdge_is_concrete():
    generic_assert_concrete_function(is_BreedingEdge)

def test_check_is_BreedingEdge_is_concrete():
    generic_assert_concrete_function(check_is_BreedingEdge)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_BreedingEdge(arch):
    assert is_BreedingEdge(arch)

def test_check_is_BreedingEdge(arch):
    with not_raises(TypeError):
        check_is_BreedingEdge(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingEdge(None, "arch")
