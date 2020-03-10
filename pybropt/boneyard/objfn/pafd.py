import numpy

def pafd(rsel, geno, wcoeff, tfreq, efreq=None, scaled=False,
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
    # generate a view of the geno matrix that only contains 'rsel' rows.
    sgeno = geno[:,rsel,:]

    # calculate the number of chromosome phases present (depth * rows)
    phases = sgeno.shape[0] * sgeno.shape[1]

    # increase precision to 'float64' for distance calculations
    if tfreq.dtype != 'float64':
        tfreq = tfreq.astype('float64')

    # calculate population frequencies in 'float64' precision
    pfreq = sgeno.sum((0,1)) / numpy.float64(phases)

    # calculate distance between target and population
    dist = numpy.absolute(tfreq - pfreq)

    # if we have scaling
    if scaled:
        # calculate edge frequencies if they haven't been provided
        if efreq is None:
            efreq = numpy.where(
                tfreq < 0.5,
                numpy.float64(1.0),
                numpy.float64(0.0)
            )
        # scale abs(target-population) by abs(target-edge)
        dist /= numpy.absolute(tfreq - efreq)

    # calculate population allele distance score using max precision
    pafd = (wcoeff * dist).sum(dtype='float64')

    # return score as the specified output data type
    return dtype.type(pafd)

# helper function
def pafd_prime(rsel, geno, wcoeff, tfreq, efreq=None,
              dtype=numpy.dtype("float64")):
    return pafd(rsel, geno, wcoeff, tfreq, efreq, scaled=True, dtype=dtype)
