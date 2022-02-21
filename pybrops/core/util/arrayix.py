"""
Module with utility functions for generating matrix indices.
"""

import numpy

__all__ = ["sqarrayix", "triuix", "triuix"]

def sqarrayix(n,k):
    """
    Generate lists containing indices for indexing square arrays.

    Parameters
    ----------
    n : int
        Length of the array along each dimension.
    k : int
        Number of dimensions of the array.

    Yields
    ------
    out : list
        List containing indices for indexing the square array.
    """
    def recurse(l,n,k):
        if len(l) == k-1:
            for i in range(n):
                l.append(i)
                yield list(l)
                l.pop()
        else:
            for i in range(n):
                l.append(i)
                yield from recurse(l,n,k)
                l.pop()
    yield from recurse([],n,k)

def triuix(n,k):
    """
    Generate lists containing indices for indexing upper triangle arrays
    including elements along the diagonal.

    Parameters
    ----------
    n : int
        Length of the array along each dimension.
    k : int
        Number of dimensions of the array.

    Yields
    ------
    out : list
        List containing indices for indexing the square array.
    """
    def recurse(l,n,k):
        st = l[-1] if len(l) else 0
        if len(l) == k-1:
            for i in range(st,n):
                l.append(i)
                yield list(l)
                l.pop()
        else:
            for i in range(st,n):
                l.append(i)
                yield from recurse(l,n,k)
                l.pop()
    yield from recurse([],n,k)

def triudix(n,k):
    """
    Generate lists containing indices for indexing upper triangle arrays
    excluding elements along the diagonal.

    Parameters
    ----------
    n : int
        Length of the array along each dimension.
    k : int
        Number of dimensions of the array.

    Yields
    ------
    out : list
        List containing indices for indexing the square array.
    """
    def recurse(l,n,k):
        st = l[-1]+1 if len(l) else 0
        if len(l) == k-1:
            for i in range(st,n):
                l.append(i)
                yield list(l)
                l.pop()
        else:
            for i in range(st,n):
                l.append(i)
                yield from recurse(l,n,k)
                l.pop()
    yield from recurse([],n,k)
