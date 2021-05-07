import inspect
import pytest

from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import KosambiMapFunction
from pybropt.popgen.gmap import is_KosambiMapFunction
from pybropt.popgen.gmap import check_is_KosambiMapFunction
from pybropt.popgen.gmap import cond_check_is_KosambiMapFunction

@pytest.fixture
def gmap(shared_datadir):
    yield ExtendedGeneticMap.from_egmap(shared_datadir / "sample.egmap")

@pytest.fixture
def kosambi():
    yield KosambiMapFunction()

def test_is_KosambiMapFunction(kosambi):
    assert is_KosambiMapFunction(kosambi)

def test_check_is_KosambiMapFunction():
    with pytest.raises(TypeError):
        check_is_KosambiMapFunction(None, "None")

def test_cond_check_is_KosambiMapFunction():
    with pytest.raises(TypeError):
        cond_check_is_KosambiMapFunction(0, "0")
