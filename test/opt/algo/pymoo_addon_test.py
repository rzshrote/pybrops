import numpy
import pytest

from pybrops.opt.algo.pymoo_addon import ReducedExchangeCrossover
# from .common_fixtures import *

@pytest.fixture
def ndecn():
    yield 20

@pytest.fixture
def nparent():
    yield 2

@pytest.fixture
def nindiv():
    yield 10

@pytest.fixture
def noption(ndecn):
    yield ndecn * 3

@pytest.fixture
def Xmat(ndecn, nparent, nindiv, noption):
    yield numpy.random.randint(0, noption, (nparent, nindiv, ndecn))

@pytest.fixture
def rex():
    yield ReducedExchangeCrossover()

### Reduced Exchange Crossover tests ###
def test_ReducedExchangeCrossover_do(rex, Xmat):
    Pmat = rex._do(None, Xmat)
    assert numpy.all(numpy.any(Xmat != Pmat, axis = 2))
