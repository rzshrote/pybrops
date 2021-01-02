import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.popgen.cmat import CoancestryMatrix

@pytest.fixture
def cmat():
    yield CoancestryMatrix()

@pytest.fixture
def vmethods(cmat):
    yield [m for m in dir(cmat) if m.startswith('__') is False]

def test_abstract_methods(cmat, vmethods):
    generic_test_abstract_methods(cmat, vmethods)
