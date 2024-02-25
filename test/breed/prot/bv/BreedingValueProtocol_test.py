import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
from pybrops.breed.prot.bv.BreedingValueProtocol import check_is_BreedingValueProtocol
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def bvprot():
    yield DummyBreedingValueProtocol()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(BreedingValueProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_method_isconcrete(BreedingValueProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_estimate_is_abstract():
    assert_method_isabstract(BreedingValueProtocol, "estimate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingValueProtocol_is_concrete():
    assert_function_isconcrete(check_is_BreedingValueProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingValueProtocol(bvprot):
    with not_raises(TypeError):
        check_is_BreedingValueProtocol(bvprot, "bvprot")
    with pytest.raises(TypeError):
        check_is_BreedingValueProtocol(None, "bvprot")
