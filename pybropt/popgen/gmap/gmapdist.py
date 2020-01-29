import numpy

def gmapdist(rgmap, cgmap):
    """
    Calculates distances between genetic map distances in 'rgmap' and 'cgmap'.
    This is done for each permutation of elements in 'rgmap' and 'cgmap'.

    Parameters
    ----------
    rgmap : numpy.ndarray
        An array of genetic map distances in cM or other linear measures of
        genetic distance representing the downstream map distances. This should
        be a 1D array. The number of elements in this array determines the
        number of rows in the output matrix.
    cgmap : numpy.ndarray
        An array of genetic map distances in cM or other linear measures of
        genetic distance representing the upstream map distances. This should be
        a 1D array. The number of elements in this array determines the number
        of columns in the output matrix.

    Returns
    -------
    d : numpy.ndarray
        A numpy.ndarray of genetic map distances. The number of rows is
        equivalent to the length of 'rgmap'. The number of columns is equivalent
        to the length of 'cgmap'.
    """
    # make an mesh grid for row and column gmap's
    mesh = numpy.meshgrid(rgmap, cgmap, indexing='ij')

    # take the difference of the two components in the grid and return abs()
    d = numpy.absolute(mesh[0] - mesh[1])

    # return distances
    return d
