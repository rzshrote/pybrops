import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.arch.GermplasmBank import GermplasmBank
from pybropt.breed.arch.GermplasmBank import is_GermplasmBank
from pybropt.breed.arch.GermplasmBank import check_is_GermplasmBank
from pybropt.breed.arch.GermplasmBank import cond_check_is_GermplasmBank

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield GermplasmBank()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(GermplasmBank)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(GermplasmBank, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(GermplasmBank, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_GermplasmBank_is_concrete():
    generic_assert_concrete_function(is_GermplasmBank)

def test_check_is_GermplasmBank_is_concrete():
    generic_assert_concrete_function(check_is_GermplasmBank)

def test_cond_check_is_GermplasmBank_is_concrete():
    generic_assert_concrete_function(cond_check_is_GermplasmBank)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_GermplasmBank(arch):
    assert is_GermplasmBank(arch)

def test_check_is_GermplasmBank(arch):
    with not_raises(TypeError):
        check_is_GermplasmBank(arch, "arch")
    with pytest.raises(TypeError):
        check_is_GermplasmBank(None, "arch")
