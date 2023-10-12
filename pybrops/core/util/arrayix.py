"""
Module with utility functions for generating matrix indices.
"""

__all__ = [
    "sqarrayix", 
    "triuix", 
    "triuix",
    "xmapix",
    "sliceaxisix",
]

import numpy


def sqarrayix(n: int, k: int) -> list:
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

def triuix(n: int, k: int) -> list:
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

def triudix(n: int, k: int) -> list:
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

def xmapix(ntaxa: int, nparent: int, unique_parents: bool) -> list:
    """
    Generate lists containing indices for cross maps.

    Parameters
    ----------
    ntaxa : int
        Number of parental candidates.
    nparent : int
        Number of parents to select from the parental candidates per cross.
    unique_parents : bool
        Whether parents should be unique. If ``False`` allow self crosses.

    Returns
    -------
    out : list
        A list of length ``nparent`` containing parental indices.
    """
    if unique_parents:
        yield from triudix(ntaxa,nparent)
    else:
        yield from triuix(ntaxa,nparent)

def sliceaxisix(shape: tuple, axis: tuple) -> tuple:
    """
    Generate tuples containing indices for iteratively slicing 
    a ``numpy.ndarray`` along a set of axes.

    Parameters
    ----------
    shape : tuple
        The shape of the ``numpy.ndarray``.
    axis : tuple
        The axes along which to iterate.
    
    Yields
    ------
    out : tuple
        A tuple which can be used to slice an array.
    """
    def recurse(l: list, s: tuple, a: tuple) -> tuple:
        if len(l) == len(s)-1:
            if len(l) in a:
                for i in range(s[len(l)]):
                    l.append(i)
                    yield tuple(l)
                    l.pop()
            else:
                l.append(slice(None))
                yield tuple(l)
                l.pop()
        else:
            if len(l) in a:
                for i in range(s[len(l)]):
                    l.append(i)
                    yield from recurse(l, s, a)
                    l.pop()
            else:
                l.append(slice(None))
                yield from recurse(l, s, a)
                l.pop()
    yield from recurse([],shape,axis)

def flattenix(
        arr: numpy.ndarray
    ) -> tuple[numpy.ndarray,tuple[numpy.ndarray,...]]:
    """
    Flatten an array and return a tuple of axis indices corresponding to their 
    position along that axis.

    Parameters
    ----------
    arr : numpy.ndarray
        An array to flatten and generate an axis index tuple.
    
    Returns
    -------
    out : tuple
        A tuple of the flattened array and axis indices.
        ``(flatarr, axisix)``
    """
    # generate a tuple of axis indices
    xi = tuple(numpy.arange(e) for e in arr.shape)

    # generate index meshgrid
    ximesh = numpy.meshgrid(*xi, indexing='ij')

    # flatten array
    flatarr = arr.flatten('C')

    # flatten index meshgrid
    axisix = tuple(mesh.flatten('C') for mesh in ximesh)

    # construct output tuple
    out = (flatarr, axisix)

    return out
