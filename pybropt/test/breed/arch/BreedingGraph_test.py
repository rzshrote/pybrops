import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import BreedingGraph
from pybropt.breed.arch import is_BreedingGraph
from pybropt.breed.arch import check_is_BreedingGraph
from pybropt.breed.arch import cond_check_is_BreedingGraph

@pytest.fixture
def bgraph():
    yield BreedingGraph()

@pytest.fixture
def vmethods(bgraph):
    yield [m for m in dir(bgraph) if m.startswith('__') is False]

def test_abstract_methods(bgraph, vmethods):
    generic_test_abstract_methods(bgraph, vmethods)

def test_is_BreedingGraph(bgraph):
    assert is_BreedingGraph(bgraph)

def test_check_is_BreedingGraph():
    with pytest.raises(TypeError):
        check_is_BreedingGraph(None, "None")

def test_cond_check_is_BreedingGraph():
    with pytest.raises(TypeError):
        cond_check_is_BreedingGraph(0, "0")
