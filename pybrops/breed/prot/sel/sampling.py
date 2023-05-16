"""
Module containing utility functions for sampling individuals.
"""

from typing import Union
import numpy
from pybrops.core.random.prng import global_prng

__all__ = [
    "stochastic_universal_sampling",
    "two_way_outcross_shuffle"
]

def stochastic_universal_sampling(
        k: int, 
        contrib: numpy.ndarray, 
        rng: Union[numpy.random.Generator,numpy.random.RandomState] = None
    ) -> numpy.ndarray:
    """
    Perform stochastic universal sampling.
    
    Parameters
    ----------
    k : int
        Number of individuals to sample.
    contrib : numpy.ndarray
        Contribution matrix of shape ``(n,)``.

        Where:

        - ``n`` is the number of individuals.

        Restrictions:

        - Values are restricted to :math:`[0,\\infty]`.
        - Sum of values in vector must be :math:`>0.0`.
    rng: Any
        Random number source. If None, use default rng.
    
    Returns
    -------
    out : numpy.ndarray
        Sampling of indices for ``contrib``.
    """
    if rng is None:                                 # if no rng provided
        rng = global_prng                           # set to global rng
    tot_fit = contrib.sum()                         # calculate the total fitness
    ptr_dist = tot_fit / k                          # calculate the distance between pointers
    indices = contrib.argsort()[::-1]               # get indices for sorted individuals
    cumsum = contrib[indices].cumsum()              # get cumulative sum of elements
    offset = rng.uniform(0.0, ptr_dist)             # get random start point
    sel = []                                        # declare output list
    ix = 0                                          # index for cumsum
    ptrs = numpy.arange(offset, tot_fit, ptr_dist)  # create pointers
    for ptr in ptrs:                                # for each pointer
        while cumsum[ix] < ptr:                     # advance index to correct location
            ix += 1
        sel.append(indices[ix])                     # append index with element
    sel = numpy.array(sel)                          # convert to ndarray
    return sel

def two_way_outcross_shuffle(
        sel: numpy.ndarray, 
        rng: Union[numpy.random.Generator,numpy.random.RandomState] = None
    ) -> numpy.ndarray:
    """
    Shuffle individuals ensuring they do not mate with themselves.

    Parameters
    ----------
    sel : numpy.ndarray
        Array of indices of individuals to select.
    rng: Any
        Random number source. If None, use default rng.
    """
    if rng is None:                                 # if no rng provided
        rng = global_prng                           # set to global rng
    rng.shuffle(sel)                                # start with random shuffle
    for st in range(0, len(sel), 2):                # for each female
        sp = st+1                                   # get stop index (male)
        j = st+2                                    # get one beyond stop index
        while (sel[st]==sel[sp]) and (j<len(sel)):  # if female == male && within array bounds
            sel[sp], sel[j] = sel[j], sel[sp]       # exchange with next entry
            j += 1                                  # increment beyond index
    return sel

