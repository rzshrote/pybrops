import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import BreedingNode
from pybropt.breed.arch import is_BreedingNode
from pybropt.breed.arch import check_is_BreedingNode
from pybropt.breed.arch import cond_check_is_BreedingNode

@pytest.fixture
def bnode():
    yield BreedingNode()

@pytest.fixture
def vmethods(bnode):
    yield [m for m in dir(bnode) if m.startswith('__') is False]

def test_abstract_methods(bnode, vmethods):
    generic_test_abstract_methods(bnode, vmethods)

def test_is_BreedingNode(bnode):
    assert is_BreedingNode(bnode)

def test_check_is_BreedingNode():
    with pytest.raises(TypeError):
        check_is_BreedingNode(None, "None")

def test_cond_check_is_BreedingNode():
    with pytest.raises(TypeError):
        cond_check_is_BreedingNode(0, "0")
