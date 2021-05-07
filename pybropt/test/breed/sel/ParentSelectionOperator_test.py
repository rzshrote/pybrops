import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.sel import ParentSelectionOperator
from pybropt.breed.sel import is_ParentSelectionOperator
from pybropt.breed.sel import check_is_ParentSelectionOperator
from pybropt.breed.sel import cond_check_is_ParentSelectionOperator

@pytest.fixture
def bedge():
    yield ParentSelectionOperator()

@pytest.fixture
def vmethods(bedge):
    yield [m for m in dir(bedge) if m.startswith('__') is False]

def test_abstract_methods(bedge, vmethods):
    generic_test_abstract_methods(bedge, vmethods)

def test_is_ParentSelectionOperator(bedge):
    assert is_ParentSelectionOperator(bedge)

def test_check_is_ParentSelectionOperator():
    with pytest.raises(TypeError):
        check_is_ParentSelectionOperator(None, "None")

def test_cond_check_is_ParentSelectionOperator():
    with pytest.raises(TypeError):
        cond_check_is_ParentSelectionOperator(0, "0")
