import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.model.gmod import GenomicModel
from pybropt.model.gmod import is_GenomicModel
from pybropt.model.gmod import check_is_GenomicModel
from pybropt.model.gmod import cond_check_is_GenomicModel

@pytest.fixture
def gmod():
    yield GenomicModel()

@pytest.fixture
def vmethods(gmod):
    yield [m for m in dir(gmod) if m.startswith('__') is False]

def test_abstract_methods(gmod, vmethods):
    generic_test_abstract_methods(gmod, vmethods)

def test_is_GenomicModel(gmod):
    assert is_GenomicModel(gmod)

def test_check_is_GenomicModel():
    with pytest.raises(TypeError):
        check_is_GenomicModel(None, "None")

def test_cond_check_is_GenomicModel():
    with pytest.raises(TypeError):
        cond_check_is_GenomicModel(0, "0")
