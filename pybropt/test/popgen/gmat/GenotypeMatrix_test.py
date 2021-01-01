import inspect
import pytest

from pybropt.popgen.gmat import GenotypeMatrix

@pytest.fixture
def gmat():
    yield GenotypeMatrix()

@pytest.fixture
def vmethods(gmat):
    yield [m for m in dir(gmat) if m.startswith('__') is False]

def test_abstract_methods(gmat, vmethods):
    for m in vmethods:
        with pytest.raises(NotImplementedError):
            fn = getattr(gmat, m)
            p = inspect.signature(fn).parameters
            kwargs = dict.fromkeys(list(p), [None] * len(p))
            fn(**kwargs)
