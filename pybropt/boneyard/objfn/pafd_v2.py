import numpy

def pafd_v2(rsel, geno, wcoeff, tfreq, efreq=None,
            dtype=numpy.dtype("float64")):
    """
    Population Allele Frequency Distance (PAFD) objective function.

    The goal is to minimize this function. Lower is better.

    This is a bare bones function. Minimal error checking.

    Parameters
    ==========
    rsel : numpy.ndarray, list, slice
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rsel' represents a single individual's row.
        If 'rsel' == slice(None, None, None) use all individuals.
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
    wcoeff : numpy.ndarray
        An array of coefficients for allele weights. Values in 'wcoeff' have
        several assumptions:
            Critical: All values should be non-negative.
            Critical: The dtype of 'wcoeff' should be a floating point.
         This array should be of shape (column,) = (N,) where 'N' represents
         number of markers.
    tfreq : None, floating, numpy.ndarray
        A floating point or numpy.ndarray of floating points representing target
        allele frequencies.
        Example:
            tfreq = numpy.array([0.2, 0.6, 0.7])
    return
    """
    # increase precision to 'float64' for distance calculations
    # use numpy.float64(a) instead of a.astype('float64') for speed
    if wcoeff.dtype != 'float64':           # if 'wcoeff' is not float64
        wcoeff = numpy.float64(wcoeff)      # cast to float64
    if tfreq.dtype != 'float64':            # if 'tfreq' is not float64
        tfreq = numpy.float64(tfreq)        # cast to float64

    # calculate the number of chromosome phases present (depth * rows)
    phases = geno[:,rsel,:].shape[0] * geno[:,rsel,:].shape[1]

    # calculate population frequencies in 'float64' precision
    pfreq = geno[:,rsel,:].sum((0,1)) / numpy.float64(phases)

    # calculate dot product of the distance (PAFD)
    pafd = numpy.absolute(tfreq - pfreq).dot(wcoeff)

    # return score as the specified output data type
    return dtype.type(pafd)
