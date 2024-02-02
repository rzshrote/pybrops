"""
Module containing iterator related tools
"""

from typing import Iterator
from itertools import tee

def iterlen(i: Iterator) -> int:
    """
    Get the length of an iterator without destroying the iterator.

    Parameters
    ----------
    i : Iterator
        Iterator for which to get the length.
    
    Returns
    -------
    out : int
        Length of the itertor.
    """
    original, copy = tee(i)
    l = sum(1 for _ in copy)
    return l
