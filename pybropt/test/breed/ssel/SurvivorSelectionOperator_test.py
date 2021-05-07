import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.ssel import SurvivorSelectionOperator
from pybropt.breed.ssel import is_SurvivorSelectionOperator
from pybropt.breed.ssel import check_is_SurvivorSelectionOperator
from pybropt.breed.ssel import cond_check_is_SurvivorSelectionOperator

@pytest.fixture
def bedge():
    yield SurvivorSelectionOperator()

@pytest.fixture
def vmethods(bedge):
    yield [m for m in dir(bedge) if m.startswith('__') is False]

def test_abstract_methods(bedge, vmethods):
    generic_test_abstract_methods(bedge, vmethods)

def test_is_SurvivorSelectionOperator(bedge):
    assert is_SurvivorSelectionOperator(bedge)

def test_check_is_SurvivorSelectionOperator():
    with pytest.raises(TypeError):
        check_is_SurvivorSelectionOperator(None, "None")

def test_cond_check_is_SurvivorSelectionOperator():
    with pytest.raises(TypeError):
        cond_check_is_SurvivorSelectionOperator(0, "0")
