import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import BreedingProgram
from pybropt.breed.arch import is_BreedingProgram
from pybropt.breed.arch import check_is_BreedingProgram
from pybropt.breed.arch import cond_check_is_BreedingProgram

@pytest.fixture
def bprog():
    yield BreedingProgram()

@pytest.fixture
def vmethods(bprog):
    yield [m for m in dir(bprog) if m.startswith('__') is False]

def test_abstract_methods(bprog, vmethods):
    generic_test_abstract_methods(bprog, vmethods)

def test_is_BreedingProgram(bprog):
    assert is_BreedingProgram(bprog)

def test_check_is_BreedingProgram():
    with pytest.raises(TypeError):
        check_is_BreedingProgram(None, "None")

def test_cond_check_is_BreedingProgram():
    with pytest.raises(TypeError):
        cond_check_is_BreedingProgram(0, "0")
