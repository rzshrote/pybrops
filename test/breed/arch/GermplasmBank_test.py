import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.arch.GermplasmBank import GermplasmBank
from pybrops.breed.arch.GermplasmBank import check_is_GermplasmBank
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyGermplasmBank()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(GermplasmBank)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(GermplasmBank, "__init__")

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
def test_check_is_GermplasmBank_is_concrete():
    assert_concrete_function(check_is_GermplasmBank)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GermplasmBank(arch):
    with not_raises(TypeError):
        check_is_GermplasmBank(arch, "arch")
    with pytest.raises(TypeError):
        check_is_GermplasmBank(None, "arch")