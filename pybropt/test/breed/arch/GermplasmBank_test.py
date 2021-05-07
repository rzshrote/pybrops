import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.breed.arch import GermplasmBank
from pybropt.breed.arch import is_GermplasmBank
from pybropt.breed.arch import check_is_GermplasmBank
from pybropt.breed.arch import cond_check_is_GermplasmBank

@pytest.fixture
def gbank():
    yield GermplasmBank()

@pytest.fixture
def vmethods(gbank):
    yield [m for m in dir(gbank) if m.startswith('__') is False]

def test_abstract_methods(gbank, vmethods):
    generic_test_abstract_methods(gbank, vmethods)

def test_is_GermplasmBank(gbank):
    assert is_GermplasmBank(gbank)

def test_check_is_GermplasmBank():
    with pytest.raises(TypeError):
        check_is_GermplasmBank(None, "None")

def test_cond_check_is_GermplasmBank():
    with pytest.raises(TypeError):
        cond_check_is_GermplasmBank(0, "0")
