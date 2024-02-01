import pytest

from pybrops.test.assert_python import not_raises
from pybrops.test.assert_python import assert_docstring
from pybrops.test.assert_python import assert_concrete_method
from pybrops.test.assert_python import assert_concrete_function

from pybrops.breed.op.init.InitializationOperator import InitializationOperator
from pybrops.breed.op.psel.ParentSelectionOperator import ParentSelectionOperator
from pybrops.breed.op.mate.MatingOperator import MatingOperator
from pybrops.breed.op.eval.EvaluationOperator import EvaluationOperator
from pybrops.breed.op.ssel.SurvivorSelectionOperator import SurvivorSelectionOperator

from pybrops.breed.arch.RecurrentSelectionBreedingProgram import RecurrentSelectionBreedingProgram
from pybrops.breed.arch.RecurrentSelectionBreedingProgram import check_is_RecurrentSelectionBreedingProgram
from .common_fixtures import *

################################################################################
################################ Test fixtures #################################
################################################################################
@pytest.fixture
def initop():
    yield DummyInitializationOperator()

@pytest.fixture
def pselop():
    yield DummyParentSelectionOperator()

@pytest.fixture
def mateop():
    yield DummyMatingOperator()

@pytest.fixture
def evalop():
    yield DummyEvaluationOperator()

@pytest.fixture
def sselop():
    yield DummySurvivorSelectionOperator()

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
    assert_docstring(RecurrentSelectionBreedingProgram)

################################################################################
############################# Test concrete methods ############################
################################################################################
def test_init_is_concrete():
    assert_concrete_method(RecurrentSelectionBreedingProgram, "__init__")

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
def test_check_is_RecurrentSelectionBreedingProgram_is_concrete():
    assert_concrete_function(check_is_RecurrentSelectionBreedingProgram)

################################################################################
######################### Test class utility functions #########################
################################################################################
def test_check_is_RecurrentSelectionBreedingProgram(arch):
    with not_raises(TypeError):
        check_is_RecurrentSelectionBreedingProgram(arch, "arch")
    with pytest.raises(TypeError):
        check_is_RecurrentSelectionBreedingProgram(None, "arch")
