import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.mate import MatingOperator
from pybropt.breed.mate import is_MatingOperator
from pybropt.breed.mate import check_is_MatingOperator
from pybropt.breed.mate import cond_check_is_MatingOperator

@pytest.fixture
def matop():
    yield MatingOperator()

@pytest.fixture
def vmethods(matop):
    yield [m for m in dir(matop) if m.startswith('__') is False]

def test_abstract_methods(matop, vmethods):
    generic_test_abstract_methods(matop, vmethods)

def test_is_MatingOperator(matop):
    assert is_MatingOperator(matop)

def test_check_is_MatingOperator():
    with pytest.raises(TypeError):
        check_is_MatingOperator(None, "None")

def test_cond_check_is_MatingOperator():
    with pytest.raises(TypeError):
        cond_check_is_MatingOperator(0, "0")
