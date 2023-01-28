import inspect
import pytest

from pybrops.test import generic_test_abstract_methods

from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol
from pybrops.breed.prot.sel.SelectionProtocol import is_SelectionProtocol
from pybrops.breed.prot.sel.SelectionProtocol import check_is_SelectionProtocol

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
