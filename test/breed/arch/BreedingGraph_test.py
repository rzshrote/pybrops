import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.arch.BreedingGraph import BreedingGraph
from pybrops.breed.arch.BreedingGraph import check_is_BreedingGraph
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyBreedingGraph()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(BreedingGraph)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(BreedingGraph, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(BreedingGraph, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingGraph_is_concrete():
    assert_concrete_function(check_is_BreedingGraph)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingGraph(arch):
    with not_raises(TypeError):
        check_is_BreedingGraph(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingGraph(None, "arch")