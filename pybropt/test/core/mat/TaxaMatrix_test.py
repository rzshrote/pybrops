import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.core.mat import TaxaMatrix
from pybropt.core.mat import is_TaxaMatrix

@pytest.fixture
def mat():
    yield TaxaMatrix()

@pytest.fixture
def vmethods(mat):
    yield [m for m in dir(mat) if m.startswith('__') is False]

def test_abstract_methods(mat, vmethods):
    generic_test_abstract_methods(mat, vmethods)

def test_is_TaxaMatrix(mat):
    assert is_TaxaMatrix(mat)
