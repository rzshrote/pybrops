import pytest
from pybrops.breed.prot.sel.SelectionProtocol import SelectionProtocol, check_is_SelectionProtocol

from pybrops.test.assert_python import assert_abstract_method

@pytest.fixture
def bedge():
    yield SelectionProtocol()

@pytest.fixture
def vmethods(bedge):
    yield [m for m in dir(bedge) if m.startswith('__') is False]

def test_abstract_methods(bedge, vmethods):
    for met in vmethods:
        assert_abstract_method(bedge, met)

def test_is_SelectionProtocol(bedge):
    assert isinstance(bedge, SelectionProtocol)

def test_check_is_SelectionProtocol():
    with pytest.raises(TypeError):
        check_is_SelectionProtocol(None, "None")


from pybrops.test.fixtures.breed.prot.sel.common_SelectionProtocol import *

def test_common(common_ncross, common_nparent):
    assert common_ncross == 10
    assert common_nparent == 2