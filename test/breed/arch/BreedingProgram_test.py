import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_class_documentation
from pybrops.test.assert_python import assert_method_isabstract
from pybrops.test.assert_python import assert_property_isabstract
from pybrops.test.assert_python import assert_method_isconcrete
from pybrops.test.assert_python import assert_function_isconcrete

from pybrops.breed.arch.BreedingProgram import BreedingProgram
from pybrops.breed.arch.BreedingProgram import check_is_BreedingProgram
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield DummyBreedingProgram()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_class_documentation(BreedingProgram)

################################################################################
############################# Test concrete methods ############################
################################################################################

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_start_genome_is_abstract():
    assert_property_isabstract(BreedingProgram, "start_genome")

def test_start_geno_is_abstract():
    assert_property_isabstract(BreedingProgram, "start_geno")

def test_start_pheno_is_abstract():
    assert_property_isabstract(BreedingProgram, "start_pheno")

def test_start_bval_is_abstract():
    assert_property_isabstract(BreedingProgram, "start_bval")

def test_start_gmod_is_abstract():
    assert_property_isabstract(BreedingProgram, "start_gmod")

def test_initop_is_abstract():
    assert_property_isabstract(BreedingProgram, "initop")

def test_pselop_is_abstract():
    assert_property_isabstract(BreedingProgram, "pselop")

def test_mateop_is_abstract():
    assert_property_isabstract(BreedingProgram, "mateop")

def test_evalop_is_abstract():
    assert_property_isabstract(BreedingProgram, "evalop")

def test_sselop_is_abstract():
    assert_property_isabstract(BreedingProgram, "sselop")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_initialize_is_abstract():
    assert_method_isabstract(BreedingProgram, "initialize")

def test_is_initialized_is_abstract():
    assert_method_isabstract(BreedingProgram, "is_initialized")

def test_reset_is_abstract():
    assert_method_isabstract(BreedingProgram, "reset")

def test_advance_is_abstract():
    assert_method_isabstract(BreedingProgram, "advance")

def test_evolve_is_abstract():
    assert_method_isabstract(BreedingProgram, "evolve")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_check_is_BreedingProgram_is_concrete():
    assert_function_isconcrete(check_is_BreedingProgram)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_BreedingProgram(arch):
    with not_raises(TypeError):
        check_is_BreedingProgram(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingProgram(None, "arch")
