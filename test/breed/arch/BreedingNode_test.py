import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.arch.BreedingNode import BreedingNode
from pybrops.breed.arch.BreedingNode import check_is_BreedingNode
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyBreedingNode()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(BreedingNode)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(BreedingNode, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_genome_is_abstract():
    assert_abstract_property(BreedingNode, "genome")

def test_geno_is_abstract():
    assert_abstract_property(BreedingNode, "geno")

def test_pheno_is_abstract():
    assert_abstract_property(BreedingNode, "pheno")

def test_bval_is_abstract():
    assert_abstract_property(BreedingNode, "bval")

def test_gmod_is_abstract():
    assert_abstract_property(BreedingNode, "gmod")

def test_t_cur_is_abstract():
    assert_abstract_property(BreedingNode, "t_cur")

def test_t_max_is_abstract():
    assert_abstract_property(BreedingNode, "t_max")


################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(BreedingNode, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingNode_is_concrete():
    assert_concrete_function(check_is_BreedingNode)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingNode(arch):
    with not_raises(TypeError):
        check_is_BreedingNode(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingNode(None, "arch")