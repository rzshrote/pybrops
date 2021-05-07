import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.model.gmod import LinearGenomicModel
from pybropt.model.gmod import is_LinearGenomicModel
from pybropt.model.gmod import check_is_LinearGenomicModel
from pybropt.model.gmod import cond_check_is_LinearGenomicModel

@pytest.fixture
def lgmod():
    yield LinearGenomicModel()

@pytest.fixture
def vmethods(lgmod):
    yield [m for m in dir(lgmod) if m.startswith('__') is False]

def test_abstract_methods(lgmod, vmethods):
    generic_test_abstract_methods(lgmod, vmethods)

def test_is_LinearGenomicModel(lgmod):
    assert is_LinearGenomicModel(lgmod)

def test_check_is_LinearGenomicModel():
    with pytest.raises(TypeError):
        check_is_LinearGenomicModel(None, "None")

def test_cond_check_is_LinearGenomicModel():
    with pytest.raises(TypeError):
        cond_check_is_LinearGenomicModel(0, "0")
