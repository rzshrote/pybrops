import inspect
import pytest

from pybropt.popgen.gmap import GeneticMap

@pytest.fixture
def gmap():
    yield GeneticMap()

@pytest.fixture
def vmethods(gmap):
    yield [m for m in dir(gmap) if m.startswith('__') is False]

def test_abstract_methods(gmap, vmethods):
    for m in vmethods:
        with pytest.raises(NotImplementedError):
            fn = getattr(gmap, m)
            p = inspect.signature(fn).parameters
            kwargs = dict.fromkeys(list(p), [None] * len(p))
            fn(**kwargs)
