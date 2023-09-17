import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.arch.BreedingEdge import BreedingEdge
from pybrops.breed.arch.BreedingEdge import check_is_BreedingEdge
from pybrops.test.breed.arch.common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyBreedingEdge()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(BreedingEdge)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(BreedingEdge, "__init__")

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
def test_check_is_BreedingEdge_is_concrete():
    assert_concrete_function(check_is_BreedingEdge)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingEdge(arch):
    with not_raises(TypeError):
        check_is_BreedingEdge(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingEdge(None, "arch")
