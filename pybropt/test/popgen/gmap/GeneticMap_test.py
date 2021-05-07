import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.popgen.gmap import GeneticMap
from pybropt.popgen.gmap import is_GeneticMap
from pybropt.popgen.gmap import check_is_GeneticMap
from pybropt.popgen.gmap import cond_check_is_GeneticMap

@pytest.fixture
def gmap():
    yield GeneticMap()

@pytest.fixture
def vmethods(gmap):
    yield [m for m in dir(gmap) if m.startswith('__') is False]

def test_abstract_methods(gmap, vmethods):
    generic_test_abstract_methods(gmap, vmethods)

def test_is_GeneticMap(gmap):
    assert is_GeneticMap(gmap)

def test_check_is_GeneticMap():
    with pytest.raises(TypeError):
        check_is_GeneticMap(None, "None")

def test_cond_check_is_GeneticMap():
    with pytest.raises(TypeError):
        cond_check_is_GeneticMap(0, "0")
