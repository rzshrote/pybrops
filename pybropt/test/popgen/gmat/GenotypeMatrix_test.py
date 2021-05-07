import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.popgen.gmat import GenotypeMatrix
from pybropt.popgen.gmat import is_GenotypeMatrix

@pytest.fixture
def gmat():
    yield GenotypeMatrix()

@pytest.fixture
def vmethods(gmat):
    yield [m for m in dir(gmat) if m.startswith('__') is False]

def test_abstract_methods(gmat, vmethods):
    generic_test_abstract_methods(gmat, vmethods)

def test_is_GenotypeMatrix(gmat):
    assert is_GenotypeMatrix(gmat)
