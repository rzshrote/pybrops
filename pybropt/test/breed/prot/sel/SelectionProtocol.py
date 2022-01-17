import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybropt.breed.prot.sel.SelectionProtocol import is_SelectionProtocol
from pybropt.breed.prot.sel.SelectionProtocol import check_is_SelectionProtocol
from pybropt.breed.prot.sel.SelectionProtocol import cond_check_is_SelectionProtocol

@pytest.fixture
def bedge():
    yield SelectionProtocol()

@pytest.fixture
def vmethods(bedge):
    yield [m for m in dir(bedge) if m.startswith('__') is False]

def test_abstract_methods(bedge, vmethods):
    generic_test_abstract_methods(bedge, vmethods)

def test_is_SelectionProtocol(bedge):
    assert is_SelectionProtocol(bedge)

def test_check_is_SelectionProtocol():
    with pytest.raises(TypeError):
        check_is_SelectionProtocol(None, "None")

def test_cond_check_is_SelectionProtocol():
    with pytest.raises(TypeError):
        cond_check_is_SelectionProtocol(0, "0")
