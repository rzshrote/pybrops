import numpy
import pytest

def generic_test_ndarray_operator(op, v, w):
    """
    Test:
        numpy.all(
            op(*v) == op(*w)
        )

    Parameters
    ----------
    op : callable
    v : tuple
    w : tuple
    """
    assert numpy.all(op(*v) == op(*w))

def generic_test_ndarray_function_eq(vfn, vargs, vkwargs, wfn, wargs, wkwargs):
    """
    Test:
        numpy.all(
            vfn(*vargs, **vkwargs) == wfn(*wargs, **wkwargs)
        )

    Parameters
    ----------
    vfn : callable
    vargs : tuple
    vkwargs : dict
    wfn : callable
    wargs : tuple
    wkwargs : dict
    """
    assert numpy.all(vfn(*vargs, **vkwargs) == wfn(*wargs, **wkwargs))
