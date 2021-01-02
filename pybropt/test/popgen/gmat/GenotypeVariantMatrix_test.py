import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.popgen.gmat import GenotypeVariantMatrix

@pytest.fixture
def gvmat():
    yield GenotypeVariantMatrix()

@pytest.fixture
def vmethods(gvmat):
    yield [m for m in dir(gvmat) if m.startswith('__') is False]

def test_abstract_methods(gvmat, vmethods):
    generic_test_abstract_methods(gvmat, vmethods)
