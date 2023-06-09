import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol
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
    assert_docstring(BreedingValueProtocol)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(BreedingValueProtocol, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_estimate_is_abstract():
    assert_abstract_method(BreedingValueProtocol, "estimate")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingValueProtocol_is_concrete():
    assert_concrete_function(check_is_BreedingValueProtocol)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingValueProtocol(bvprot):
    with not_raises(TypeError):
        check_is_BreedingValueProtocol(bvprot, "bvprot")
    with pytest.raises(TypeError):
        check_is_BreedingValueProtocol(None, "bvprot")
