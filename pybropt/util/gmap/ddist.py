import numpy

def ddist(d0, d1):
    """
    Calculates distances between genetic map distances in 'd0' and 'd1'. This is
    done for each permutation of elements in 'd0' and 'd1'. The difference is
    calculated by abs(d1 - d0).

    Parameters
    ----------
    d0 : numpy.ndarray
        An array of genetic map distances in cM or other linear measures of
        genetic distance representing the upstream map distances. This should be
        a 1D array. The number of elements in this array determines the number
        of columns in the output matrix.
    d1 : numpy.ndarray
        An array of genetic map distances in cM or other linear measures of
        genetic distance representing the downstream map distances. This should
        be a 1D array. The number of elements in this array determines the
        number of rows in the output matrix.

    Returns
    -------
    d : numpy.ndarray
        Returns an array with (rows, columns) = (len(d1), len(d0)).
    """
    # TODO: test tile vs repeat speeds.
    # TODO: If if find that I keep on having to transform matrices from this
    #       function, I'll change the output shape.

    # get rows and columns of the expected output matrix
    # TODO: See if storing in r, c is faster than calling the array.
    r = d1.shape[0]
    c = d0.shape[0]

    # first start by tiling the the d0 matrix to a dimension of
    # (len(d1), len(d0))
    # NOTE: this assumes C style matrix layout.
    U = numpy.tile(d0, (r, 1))

    # next repeat the d1 array and reshape it to be the same size as U.
    # NOTE: this assumes C style matrix layout.
    D = numpy.repeat(d1, c)
    D.shape = (r, c)

    # take the absolute value of the difference between R and its transpose.
    return numpy.absolute(D - U)
