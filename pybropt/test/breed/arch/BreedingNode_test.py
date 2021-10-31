import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.arch import BreedingNode
from pybropt.breed.arch import is_BreedingNode
from pybropt.breed.arch import check_is_BreedingNode
from pybropt.breed.arch import cond_check_is_BreedingNode

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield BreedingNode()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(BreedingNode)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(BreedingNode, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(BreedingNode, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_BreedingNode_is_concrete():
    generic_assert_concrete_function(is_BreedingNode)

def test_check_is_BreedingNode_is_concrete():
    generic_assert_concrete_function(check_is_BreedingNode)

def test_cond_check_is_BreedingNode_is_concrete():
    generic_assert_concrete_function(cond_check_is_BreedingNode)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_BreedingNode(arch):
    assert is_BreedingNode(arch)

def test_check_is_BreedingNode(arch):
    with not_raises(TypeError):
        check_is_BreedingNode(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingNode(None, "arch")
