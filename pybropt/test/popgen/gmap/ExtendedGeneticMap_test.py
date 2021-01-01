import pytest
from pybropt.popgen.gmap import ExtendedGeneticMap
from pybropt.popgen.gmap import is_ExtendedGeneticMap

@pytest.fixture
def egmap(shared_datadir):
    data_path = shared_datadir / "McMullen_2009_US_NAM.M.egmap"
    yield ExtendedGeneticMap.from_egmap(data_path)

def test_type(egmap):
    assert is_ExtendedGeneticMap(egmap)
