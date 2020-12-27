import numpy
from itertools import chain

# generator to calculate pairwise products in chunks
def pweight(coeff, start, stop, divisor = 1.0, chunk = None):
    """
    Generator function for calulating pairwise marker weights scaled by
    'divisor'.

    Parameters
    ==========
    coeff : numpy.ndarray
        1D array of marker coefficients.
    start : int
        Starting index to compute matrix chunks.
    stop : int
        Stopping index to comput matrix chunks.
    divisor : float, numpy.float
        Divisor to divide the results by.
    chunk : None, int
        Number of markers to process at once.

    Yields
    ======
    yield: numpy.ndarray
        A numpy.ndarray of marker coefficients for a matrix chunk.
    """
    # if the chunk size is None, set chunk to stop-start (compute whole matrix)
    if chunk is None:
        chunk = stop - start

    # for each row, and column get row start and stops
    # gory details:
    #     # this calculates start indices: [start, stop) by chunk
    #     rst = range(start, stop, chunk)
    #     # this calculates stop indices: [start+chunk, stop] by chunk
    #     rsp = chain(
    #               # this calculates all but one index: [chunk, l) by chunk
    #               range(start+chunk, stop, chunk),
    #               # tack on the last index by making into tuple
    #               (stop,)
    #           )
    #     # finally zip everything up

    # TODO: Potential optimization where we compute numpy.triu fragments first,
    #       then the rest of the matrix. This would save an if else at the cost
    #       of additional iterators. Not sure if this is worth it.
    for rst,rsp in zip(                     # zip start, stop indices
        range(start,stop,chunk),            # calculate start indices
        chain(                              # calculate stop indices
            range(start+chunk,stop,chunk),  # calculate first several
            (stop,)                         # tack on last index
        )
    ):
        for cst,csp in zip(                     # zip start, stop indices
            range(start,stop,chunk),            # calculate start indices
            chain(                              # calculate stop indices
                range(start+chunk,stop,chunk),  # calculate first several
                (stop,)                         # tack on last index
            )
        ):
            yield (coeff[rst:rsp,None] @ coeff[None,cst:csp]) / divisor

# a generator to create an upper triangle matrix of pairwise products in chunks
def pweight_triu(coeff, start, stop, divisor = 1.0, chunk = None):
    """
    Generator function for calulating an upper triangle matrix of pairwise
    marker weights scaled by 'divisor'.

    Parameters
    ==========
    coeff : numpy.ndarray
        1D array of marker coefficients.
    start : int
        Starting index to compute matrix chunks.
    stop : int
        Stopping index to comput matrix chunks.
    divisor : float, numpy.float
        Divisor to divide the results by.
    chunk : None, int
        Number of elements in coeff to process at once.

    Yields
    ======
    yield: numpy.ndarray
        A numpy.ndarray of marker coefficients for a matrix chunk.
    """
    # if the chunk size is None, set chunk to stop-start (compute whole matrix)
    if chunk is None:
        chunk = stop - start

    # for each row, and column get row start and stops
    # gory details:
    #     # this calculates start indices: [start, stop) by chunk
    #     rst = range(start, stop, chunk)
    #     # this calculates stop indices: [start+chunk, stop] by chunk
    #     rsp = chain(
    #               # this calculates all but one index: [chunk, l) by chunk
    #               range(start+chunk, stop, chunk),
    #               # tack on the last index by making into tuple
    #               (stop,)
    #           )
    #     # finally zip everything up

    # TODO: Potential optimization where we compute numpy.triu fragments first,
    #       then the rest of the matrix. This would save an if else at the cost
    #       of additional iterators. Not sure if this is worth it.
    for rst,rsp in zip(                     # zip start, stop indices
        range(start,stop,chunk),            # calculate start indices
        chain(                              # calculate stop indices
            range(start+chunk,stop,chunk),  # calculate first several
            (stop,)                         # tack on last index
        )
    ):
        for cst,csp in zip(                     # zip start, stop indices
            range(rst,stop,chunk),              # calculate start indices
            chain(                              # calculate stop indices
                range(rst+chunk,stop,chunk),    # calculate first several
                (stop,)                         # tack on last index
            )
        ):
            # if the row start and column start are not equal, then the matrix
            # chunk is off the diagonal and we return the whole thing.
            if rst != cst:
                yield (coeff[rst:rsp,None] @ coeff[None,cst:csp]) / divisor
            else:
                yield numpy.triu(coeff[rst:rsp,None] @ coeff[None,cst:csp]) / divisor
