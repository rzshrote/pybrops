import inspect
import pytest

from pybropt.test import generic_test_abstract_methods

from pybropt.popgen.gmap import GeneticMapFunction
from pybropt.popgen.gmap import is_GeneticMapFunction
from pybropt.popgen.gmap import check_is_GeneticMapFunction
from pybropt.popgen.gmap import cond_check_is_GeneticMapFunction

@pytest.fixture
def gmapfn():
    yield GeneticMapFunction()

@pytest.fixture
def vmethods(gmapfn):
    yield [m for m in dir(gmapfn) if m.startswith('__') is False]

def test_abstract_methods(gmapfn, vmethods):
    generic_test_abstract_methods(gmapfn, vmethods)

def test_is_GeneticMapFunction(gmapfn):
    assert is_GeneticMapFunction(gmapfn)

def test_check_is_GeneticMapFunction():
    with pytest.raises(TypeError):
        check_is_GeneticMapFunction(None, "None")

def test_cond_check_is_GeneticMapFunction():
    with pytest.raises(TypeError):
        cond_check_is_GeneticMapFunction(0, "0")
