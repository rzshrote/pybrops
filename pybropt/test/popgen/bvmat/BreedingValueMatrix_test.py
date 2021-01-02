import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.popgen.bvmat import BreedingValueMatrix

@pytest.fixture
def bvmat():
    yield BreedingValueMatrix()

@pytest.fixture
def vmethods(bvmat):
    yield [m for m in dir(bvmat) if m.startswith('__') is False]

def test_abstract_methods(bvmat, vmethods):
    generic_test_abstract_methods(bvmat, vmethods)
