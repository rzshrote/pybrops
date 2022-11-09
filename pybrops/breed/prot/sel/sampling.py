"""
Module containing utility functions for sampling individuals.
"""

import numpy
from pybrops.core.random import global_prng

def stochastic_universal_sampling(k: int, contrib: numpy.ndarray, rng = None):
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
