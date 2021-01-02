import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.popgen.gmat import GenotypeMatrix

@pytest.fixture
def gmat():
    yield GenotypeMatrix()

@pytest.fixture
def vmethods(gmat):
    yield [m for m in dir(gmat) if m.startswith('__') is False]

def test_abstract_methods(gmat, vmethods):
    generic_test_abstract_methods(gmat, vmethods)
