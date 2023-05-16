import inspect
import pytest

from pybrops.test.assert_python import assert_abstract_methods

from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import UnconstrainedSelectionProtocol
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import is_SelectionProtocol
from pybrops.breed.prot.sel.UnconstrainedSelectionProtocol import check_is_SelectionProtocol

@pytest.fixture
def bedge():
    yield UnconstrainedSelectionProtocol()

@pytest.fixture
def vmethods(bedge):
    yield [m for m in dir(bedge) if m.startswith('__') is False]

def test_abstract_methods(bedge, vmethods):
    assert_abstract_methods(bedge, vmethods)

def test_is_SelectionProtocol(bedge):
    assert is_SelectionProtocol(bedge)

def test_check_is_SelectionProtocol():
    with pytest.raises(TypeError):
        check_is_SelectionProtocol(None, "None")
