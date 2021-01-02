import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.core.mat import Matrix
from pybropt.core.mat import is_Matrix

@pytest.fixture
def mat():
    yield Matrix()

@pytest.fixture
def vmethods(mat):
    yield [m for m in dir(mat) if m.startswith('__') is False]

def test_abstract_methods(mat, vmethods):
    generic_test_abstract_methods(mat, vmethods)

def test_is_Matrix(mat):
    assert is_Matrix(mat)
