"""
Module containing sampling subroutines.
"""

__all__ = [
    "stochastic_universal_sampling",
    "tiled_choice",
    "axis_shuffle",
    "outcross_shuffle"
]

from numbers import Integral
from typing import Optional, Tuple, Union
import numpy
from numpy.random import Generator, RandomState

from pybrops.core.random.prng import global_prng
from pybrops.core.util.arrayix import sliceaxisix


def stochastic_universal_sampling(
        a: numpy.ndarray,
        p: numpy.ndarray, 
        size: Optional[Union[Integral,Tuple[Integral,...]]] = None,
        rng: Optional[Union[Generator,RandomState]] = None
    ) -> numpy.ndarray:
    """
    Perform stochastic universal sampling.
    
    Parameters
    ----------
    a : numpy.ndarray
        Array of shape ``(n,)`` containing element from which to sample.
    p : numpy.ndarray
        Contribution probability array of shape ``(n,)``.

        Where:

        - ``n`` is the number of elements.

        Restrictions:

        - Values are restricted to :math:`[0,\\infty]`.
        - Sum of values in vector must be :math:`>0.0`.
    size : Integral, Tuple[Integral,...]
        Target output size.
    rng: Generator, RandomState, None
        Random number source. If None, use default rng.
    
    Returns
    -------
    out : numpy.ndarray
        Sampling of items from ``a``.
    """
    if rng is None:                                 # if no rng provided
        rng = global_prng                           # set to global rng
    if isinstance(size, Integral):                  # if size is integer
        size = (size,)                              # convert size to tuple
    k = numpy.product(size)                         # calculate total number of elements we need
    tot_fit = p.sum()                               # calculate the total fitness
    ptr_dist = tot_fit / k                          # calculate the distance between pointers
    indices = p.argsort()[::-1]                     # get indices for sorted individuals
    cumsum = p[indices].cumsum()                    # get cumulative sum of elements
    offset = rng.uniform(0.0, ptr_dist)             # get random start point
    sel = []                                        # declare output list
    ix = 0                                          # index for cumsum
    ptrs = numpy.arange(offset, tot_fit, ptr_dist)  # create pointers
    for ptr in ptrs:                                # for each pointer
        while cumsum[ix] < ptr:                     # advance index to correct location
            ix += 1
        sel.append(indices[ix])                     # append index with element
    sel = numpy.array(sel)                          # convert to ndarray
    rng.shuffle(sel)                                # shuffle the output array
    sel = sel.reshape(size)                         # reshape array
    return a[sel]

def tiled_choice(
        a: numpy.ndarray,
        size: Optional[Union[Integral,Tuple[Integral,...]]] = None,
        replace: bool = True,
        p: Optional[numpy.ndarray] = None,
        rng: Optional[Union[Generator,RandomState]] = None
    ) -> numpy.ndarray:
    """
    Generate a random sample from a given 1-D array.

    If sampling without replacement, tile the sampling by filling with complete 
    sets until a complete set cannot be put into the output, then sample without 
    replacement to fill the remaining samples in the output.

    Parameters
    ----------
    a : numpy.ndarray, Integral
        If an ndarray, a random sample is generated from its elements.
        If an Integral, the random sample is generated as if it were np.arange(a)
    size : Integral, Tuple of Integral, None
        Output shape. If the given shape is, e.g., (m, n, k), then m * n * k samples are drawn.
        Default is None, in which case a single value is returned.
    replace : bool
        Whether the sample is with or without replacement.
        Default is True, meaning that a value of a can be selected multiple times.
    p : numpy.ndarray, None
        The probabilities associated with each entry in a. 
        If not given, the sample assumes a uniform distribution over all entries in a.
    rng : numpy.random.Generator, numpy.random.RandomState, None
        A random number generator source.
        If None, then use the default global random number generator.

    Returns
    -------
    out : numpy.ndarray
        The generated random samples
    """
    # convert size to tuple
    if isinstance(size, Integral):
        size = (size,)
    
    # ensure we have a random number generators source
    if rng is None:
        rng = global_prng
    
    # get shape parameters
    nsample = numpy.product(size)

    # allocate memory for output
    out = numpy.empty(nsample, dtype = a.dtype)

    # if sampling with replacement, then tiling is not needed
    if replace:
        out = rng.choice(a, size, replace, p)
    
    # if sampling without replacement, then we need to tile
    else:
        # get shape parameters
        noption = len(a)

        # get quotient and remainder
        qu, re = divmod(nsample, noption)

        # tile the decisions
        for i in range(qu):
            out[(i*noption):((i+1)*noption)] = a
        
        # fill the remaining decisions with random choice
        out[(qu*noption):] = rng.choice(a, re, replace, p)

        # shuffle the decisions
        rng.shuffle(out)

        # reshape output
        out = out.reshape(size)
    
    return out

def axis_shuffle(
        a: numpy.ndarray,
        axis: Optional[Union[Integral,Tuple[Integral,...]]] = None,
        rng: Optional[Union[Generator,RandomState]] = None
    ) -> None:
    """
    Perform an in-place shuffle for an array along a set of axes.

    Parameters
    ----------
    a : numpy.ndarray
        An array to shuffle.
    axis : Integral, Tuple[Integral,...]
        Axis or axes along which to perform an in-place shuffle.
    rng : numpy.random.Generator, numpy.random.RandomState, None
        A random number generator source.
        If None, then use the default global random number generator.
    """
    # check type for ``x``
    if not isinstance(a, numpy.ndarray):
        raise TypeError("'x' must be of type numpy.ndarray, but received type '{0}'".format(type(a).__name__))
    
    # check type for ``axis``
    if isinstance(axis, Integral):
        axis = (axis,)
    elif isinstance(axis, tuple):
        pass
    else:
        raise TypeError("'axis' must be an Integral type or a tuple, but recieved type '{0}'".format(type(axis).__name__))
    
    # check type for ``rng``
    if rng is None:
        rng = global_prng
    if not isinstance(rng, Generator) and not isinstance(rng, RandomState):
        raise TypeError("'rng' must be of type Generator or RandomState, but received type '{0}'".format(type(rng).__name__))
    
    # main code
    for s in sliceaxisix(a.shape,axis):
        rng.shuffle(a[s])

def outcross_shuffle(
        xconfig: numpy.ndarray,
        rng: Optional[Union[Generator,RandomState]] = None
    ) -> None:
    """
    Shuffle individuals between crosses to ensure that matings are (at least 
    locally) maximally outcrossing. This function uses a stochastic descent 
    hillclimber to alternate individuals between cross configurations to 
    maximize the outcrossing (minimize identical individuals within a single 
    cross configuration).

    Parameters
    ----------
    xconfig : numpy.ndarray
        Array of shape (ncross, nparent) containing indices of individuals to shuffle.
    rng : numpy.random.Generator, numpy.random.RandomState, None
        A random number generator source.
        If None, then use the default global random number generator.
    """
    # handle prng
    if rng is None:
        rng = global_prng

    # create objective function to be minimized: the number of duplicates
    def objfn(x: numpy.ndarray):
        out = 0
        for xrow in x:
            u,c = numpy.unique(xrow, return_counts=True)
            out += numpy.sum(c-1)
        return out

    # ravel the output array to a 1D view
    xravel = xconfig.ravel()

    # get starting solution score
    gbest_score = objfn(xconfig)       # score solution

    # create exchange indices
    exchix = numpy.array([[i,j] for i in range(len(xravel)) for j in range(i+1, len(xravel))])

    # establish stopping criterion
    iterate = True

    # main loop
    while iterate:
        rng.shuffle(exchix)                     # shuffle exchange indices
        local_optima = True                     # whether we went through all exchange indices
        for i,j in exchix:                      # for each exchange pair
            xravel[i], xravel[j] = xravel[j], xravel[i] # exchange values
            score = objfn(xconfig)              # score configuration
            if score < gbest_score:             # if weighted score is better
                gbest_score = score             # store new best score
                local_optima = False            # indicate we are not in a local optima
                break                           # break from for loop
            xravel[i], xravel[j] = xravel[j], xravel[i] # exchange values back to original
        iterate = not local_optima              # update whether to continue climbing

