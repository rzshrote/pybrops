import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.model.gmod import NonlinearGenomicModel
from pybropt.model.gmod import is_NonlinearGenomicModel
from pybropt.model.gmod import check_is_NonlinearGenomicModel
from pybropt.model.gmod import cond_check_is_NonlinearGenomicModel

@pytest.fixture
def lgmod():
    yield NonlinearGenomicModel()

@pytest.fixture
def vmethods(lgmod):
    yield [m for m in dir(lgmod) if m.startswith('__') is False]

def test_abstract_methods(lgmod, vmethods):
    generic_test_abstract_methods(lgmod, vmethods)

def test_is_NonlinearGenomicModel(lgmod):
    assert is_NonlinearGenomicModel(lgmod)

def test_check_is_NonlinearGenomicModel():
    with pytest.raises(TypeError):
        check_is_NonlinearGenomicModel(None, "None")

def test_cond_check_is_NonlinearGenomicModel():
    with pytest.raises(TypeError):
        cond_check_is_NonlinearGenomicModel(0, "0")
