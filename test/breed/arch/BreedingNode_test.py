import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

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
    assert_class_documentation(BreedingNode)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_genome_is_abstract():
    assert_property_isabstract(BreedingNode, "genome")

def test_geno_is_abstract():
    assert_property_isabstract(BreedingNode, "geno")

def test_pheno_is_abstract():
    assert_property_isabstract(BreedingNode, "pheno")

def test_bval_is_abstract():
    assert_property_isabstract(BreedingNode, "bval")

def test_gmod_is_abstract():
    assert_property_isabstract(BreedingNode, "gmod")

def test_t_cur_is_abstract():
    assert_property_isabstract(BreedingNode, "t_cur")

def test_t_max_is_abstract():
    assert_property_isabstract(BreedingNode, "t_max")


################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     assert_method_isabstract(BreedingNode, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingNode_is_concrete():
    assert_function_isconcrete(check_is_BreedingNode)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingNode(arch):
    with not_raises(TypeError):
        check_is_BreedingNode(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingNode(None, "arch")
