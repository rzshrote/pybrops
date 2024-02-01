import pytest

from pybrops.core.mat.util import get_axis

def test_get_axis():
    assert get_axis(-3, 3) == 0
    assert get_axis(-2, 3) == 1
    assert get_axis(-1, 3) == 2
    assert get_axis(0, 3) == 0
    assert get_axis(1, 3) == 1
    assert get_axis(2, 3) == 2

def test_get_axis_error_low():
    with pytest.raises(IndexError):
        get_axis(-4, 3)

def test_get_axis_error_high():
    with pytest.raises(IndexError):
        get_axis(3, 3)
