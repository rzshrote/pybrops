import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.test import not_raises
from pybropt.test import generic_assert_docstring
from pybropt.test import generic_assert_abstract_method
from pybropt.test import generic_assert_abstract_function
from pybropt.test import generic_assert_abstract_property
from pybropt.test import generic_assert_concrete_method
from pybropt.test import generic_assert_concrete_function

from pybropt.breed.op.init import InitializationOperator
from pybropt.breed.op.psel import ParentSelectionOperator
from pybropt.breed.op.mate import MatingOperator
from pybropt.breed.op.eval import EvaluationOperator
from pybropt.breed.op.ssel import SurvivorSelectionOperator

from pybropt.breed.arch import RecurrentSelectionBreedingProgram
from pybropt.breed.arch import is_RecurrentSelectionBreedingProgram
from pybropt.breed.arch import check_is_RecurrentSelectionBreedingProgram
from pybropt.breed.arch import cond_check_is_RecurrentSelectionBreedingProgram

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def initop():
    yield InitializationOperator()

@pytest.fixture
def pselop():
    yield ParentSelectionOperator()

@pytest.fixture
def mateop():
    yield MatingOperator()

@pytest.fixture
def evalop():
    yield EvaluationOperator()

@pytest.fixture
def sselop():
    yield SurvivorSelectionOperator()

@pytest.fixture
def t_max():
    yield 20

@pytest.fixture
def arch(initop, pselop, mateop, evalop, sselop, t_max):
    yield RecurrentSelectionBreedingProgram(
        initop = initop,
        pselop = pselop,
        mateop = mateop,
        evalop = evalop,
        sselop = sselop,
        t_max = t_max
    )

################################################################################
############################## Test class docstring ############################
################################################################################
def test_class_docstring():
    generic_assert_docstring(RecurrentSelectionBreedingProgram)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    generic_assert_concrete_method(RecurrentSelectionBreedingProgram, "__init__")

################################################################################
########################### Test abstract properties ###########################
################################################################################

################################################################################
############################# Test abstract methods ############################
################################################################################
# def test__is_abstract():
#     generic_assert_abstract_method(RecurrentSelectionBreedingProgram, "")

################################################################################
################### Test for conrete class utility functions ###################
################################################################################
def test_is_RecurrentSelectionBreedingProgram_is_concrete():
    generic_assert_concrete_function(is_RecurrentSelectionBreedingProgram)

def test_check_is_RecurrentSelectionBreedingProgram_is_concrete():
    generic_assert_concrete_function(check_is_RecurrentSelectionBreedingProgram)

def test_cond_check_is_RecurrentSelectionBreedingProgram_is_concrete():
    generic_assert_concrete_function(cond_check_is_RecurrentSelectionBreedingProgram)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_is_RecurrentSelectionBreedingProgram(arch):
    assert is_RecurrentSelectionBreedingProgram(arch)

def test_check_is_RecurrentSelectionBreedingProgram(arch):
    with not_raises(TypeError):
        check_is_RecurrentSelectionBreedingProgram(arch, "arch")
    with pytest.raises(TypeError):
        check_is_RecurrentSelectionBreedingProgram(None, "arch")
