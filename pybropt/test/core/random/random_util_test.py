import inspect
import pytest
from numpy.random import Generator

from pybropt.core.random import spawn
from pybropt.core.random import seed

def test_spawn_None():
    s = spawn()
    assert isinstance(s, Generator)

def test_spawn_int():
    l = spawn(5)
    assert isinstance(l, list)
    for e in l:
        assert isinstance(e, Generator)

def test_seed_None():
    seed(None)

def test_seed_int():
    seed(int(123456789))
