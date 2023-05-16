import pytest

from pybrops.test.assert_python import assert_abstract_methods
from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_abstract_method
from pybrops.test.assert_python import assert_abstract_function
from pybrops.test.assert_python import assert_abstract_property
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.arch.BreedingProgram import BreedingProgram
from pybrops.breed.arch.BreedingProgram import is_BreedingProgram
from pybrops.breed.arch.BreedingProgram import check_is_BreedingProgram

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def arch():
    yield BreedingProgram()

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    assert_docstring(BreedingProgram)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(BreedingProgram, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################
def test_start_genome_is_abstract():
    assert_abstract_property(BreedingProgram, "start_genome")

def test_start_geno_is_abstract():
    assert_abstract_property(BreedingProgram, "start_geno")

def test_start_pheno_is_abstract():
    assert_abstract_property(BreedingProgram, "start_pheno")

def test_start_bval_is_abstract():
    assert_abstract_property(BreedingProgram, "start_bval")

def test_start_gmod_is_abstract():
    assert_abstract_property(BreedingProgram, "start_gmod")

def test_initop_is_abstract():
    assert_abstract_property(BreedingProgram, "initop")

def test_pselop_is_abstract():
    assert_abstract_property(BreedingProgram, "pselop")

def test_mateop_is_abstract():
    assert_abstract_property(BreedingProgram, "mateop")

def test_evalop_is_abstract():
    assert_abstract_property(BreedingProgram, "evalop")

def test_sselop_is_abstract():
    assert_abstract_property(BreedingProgram, "sselop")

################################################################################
############################# Test abstract methods ############################
################################################################################
def test_initialize_is_abstract():
    assert_abstract_method(BreedingProgram, "initialize")

def test_is_initialized_is_abstract():
    assert_abstract_method(BreedingProgram, "is_initialized")

def test_reset_is_abstract():
    assert_abstract_method(BreedingProgram, "reset")

def test_advance_is_abstract():
    assert_abstract_method(BreedingProgram, "advance")

def test_evolve_is_abstract():
    assert_abstract_method(BreedingProgram, "evolve")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_BreedingProgram_is_concrete():
    assert_concrete_function(is_BreedingProgram)

def test_check_is_BreedingProgram_is_concrete():
    assert_concrete_function(check_is_BreedingProgram)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_BreedingProgram(arch):
    assert is_BreedingProgram(arch)

def test_check_is_BreedingProgram(arch):
    with not_raises(TypeError):
        check_is_BreedingProgram(arch, "arch")
    with pytest.raises(TypeError):
        check_is_BreedingProgram(None, "arch")
