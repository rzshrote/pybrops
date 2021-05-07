import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import EmigrationOperator
from pybropt.breed.arch import is_EmigrationOperator
from pybropt.breed.arch import check_is_EmigrationOperator
from pybropt.breed.arch import cond_check_is_EmigrationOperator

@pytest.fixture
def emiop():
    yield EmigrationOperator()

@pytest.fixture
def vmethods(emiop):
    yield [m for m in dir(emiop) if m.startswith('__') is False]

def test_abstract_methods(emiop, vmethods):
    generic_test_abstract_methods(emiop, vmethods)

def test_is_EmigrationOperator(emiop):
    assert is_EmigrationOperator(emiop)

def test_check_is_EmigrationOperator():
    with pytest.raises(TypeError):
        check_is_EmigrationOperator(None, "None")

def test_cond_check_is_EmigrationOperator():
    with pytest.raises(TypeError):
        cond_check_is_EmigrationOperator(0, "0")
