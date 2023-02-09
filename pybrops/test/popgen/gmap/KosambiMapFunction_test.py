import inspect
import pytest

from pybrops.popgen.gmap.ExtendedGeneticMap import ExtendedGeneticMap
from pybrops.popgen.gmap.KosambiMapFunction import KosambiMapFunction
from pybrops.popgen.gmap.KosambiMapFunction import is_KosambiMapFunction
from pybrops.popgen.gmap.KosambiMapFunction import check_is_KosambiMapFunction

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
