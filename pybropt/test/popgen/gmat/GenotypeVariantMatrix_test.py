import inspect
import pytest

from pybropt.popgen.gmat import GenotypeVariantMatrix

@pytest.fixture
def gvmat():
    yield GenotypeVariantMatrix()

@pytest.fixture
def vmethods(gvmat):
    yield [m for m in dir(gvmat) if m.startswith('__') is False]

def test_abstract_methods(gvmat, vmethods):
    for m in vmethods:
        with pytest.raises(NotImplementedError):
            fn = getattr(gvmat, m)
            p = inspect.signature(fn).parameters
            kwargs = dict.fromkeys(list(p), [None] * len(p))
            fn(**kwargs)
