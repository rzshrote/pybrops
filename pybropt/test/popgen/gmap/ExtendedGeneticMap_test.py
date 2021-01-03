import inspect
import pytest

from pybropt.breed.arch import ExtendedGeneticMap
from pybropt.breed.arch import is_ExtendedGeneticMap
from pybropt.breed.arch import check_is_ExtendedGeneticMap
from pybropt.breed.arch import cond_check_is_ExtendedGeneticMap

@pytest.fixture
def egmap(shared_datadir):
    data_path = shared_datadir / "McMullen_2009_US_NAM.M.egmap"
    yield ExtendedGeneticMap.from_egmap(data_path)

def test_is_ExtendedGeneticMap(egmap):
    assert is_ExtendedGeneticMap(egmap)

def test_check_is_ExtendedGeneticMap():
    with pytest.raises(TypeError):
        check_is_ExtendedGeneticMap(None, "None")

def test_cond_check_is_ExtendedGeneticMap():
    with pytest.raises(TypeError):
        cond_check_is_ExtendedGeneticMap(0, "0")
