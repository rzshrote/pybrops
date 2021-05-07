import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import ImmigrationOperator
from pybropt.breed.arch import is_ImmigrationOperator
from pybropt.breed.arch import check_is_ImmigrationOperator
from pybropt.breed.arch import cond_check_is_ImmigrationOperator

@pytest.fixture
def immiop():
    yield ImmigrationOperator()

@pytest.fixture
def vmethods(immiop):
    yield [m for m in dir(immiop) if m.startswith('__') is False]

def test_abstract_methods(immiop, vmethods):
    generic_test_abstract_methods(immiop, vmethods)

def test_is_ImmigrationOperator(immiop):
    assert is_ImmigrationOperator(immiop)

def test_check_is_ImmigrationOperator():
    with pytest.raises(TypeError):
        check_is_ImmigrationOperator(None, "None")

def test_cond_check_is_ImmigrationOperator():
    with pytest.raises(TypeError):
        cond_check_is_ImmigrationOperator(0, "0")
