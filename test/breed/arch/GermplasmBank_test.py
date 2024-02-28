import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

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
    assert_class_documentation(GermplasmBank)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     assert_method_isabstract(GermplasmBank, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_GermplasmBank_is_concrete():
    assert_function_isconcrete(check_is_GermplasmBank)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_GermplasmBank(arch):
    with not_raises(TypeError):
        check_is_GermplasmBank(arch, "arch")
    with pytest.raises(TypeError):
        check_is_GermplasmBank(None, "arch")
