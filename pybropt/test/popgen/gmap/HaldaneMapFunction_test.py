import inspect
import pytest

from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import HaldaneMapFunction
from pybropt.popgen.gmap import is_HaldaneMapFunction
from pybropt.popgen.gmap import check_is_HaldaneMapFunction
from pybropt.popgen.gmap import cond_check_is_HaldaneMapFunction

@pytest.fixture
def gmap(shared_datadir):
    yield ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")

@pytest.fixture
def haldane():
    yield HaldaneMapFunction()

def test_is_HaldaneMapFunction(haldane):
    assert is_HaldaneMapFunction(haldane)

def test_check_is_HaldaneMapFunction():
    with pytest.raises(TypeError):
        check_is_HaldaneMapFunction(None, "None")

def test_cond_check_is_HaldaneMapFunction():
    with pytest.raises(TypeError):
        cond_check_is_HaldaneMapFunction(0, "0")
