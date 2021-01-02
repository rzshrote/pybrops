import inspect
import pytest

from pybropt.test import generic_test_abstract_methods
from pybropt.core.mat import SortableMatrix
from pybropt.core.mat import is_SortableMatrix

@pytest.fixture
def mat():
    yield SortableMatrix()

@pytest.fixture
def vmethods(mat):
    yield [m for m in dir(mat) if m.startswith('__') is False]

def test_abstract_methods(mat, vmethods):
    generic_test_abstract_methods(mat, vmethods)

def test_is_SortableMatrix(mat):
    assert is_SortableMatrix(mat)
