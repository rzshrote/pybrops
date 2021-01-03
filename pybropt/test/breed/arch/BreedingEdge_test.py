import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import BreedingEdge
from pybropt.breed.arch import is_BreedingEdge
from pybropt.breed.arch import check_is_BreedingEdge
from pybropt.breed.arch import cond_check_is_BreedingEdge

@pytest.fixture
def bedge():
    yield BreedingEdge()

@pytest.fixture
def vmethods(bedge):
    yield [m for m in dir(bedge) if m.startswith('__') is False]

def test_abstract_methods(bedge, vmethods):
    generic_test_abstract_methods(bedge, vmethods)

def test_is_BreedingEdge(bedge):
    assert is_BreedingEdge(bedge)

def test_check_is_BreedingEdge():
    with pytest.raises(TypeError):
        check_is_BreedingEdge(None, "None")

def test_cond_check_is_BreedingEdge():
    with pytest.raises(TypeError):
        cond_check_is_BreedingEdge(0, "0")
