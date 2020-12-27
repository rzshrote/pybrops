import numpy

def paa(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64")):
    """
    Population Allele Availability (PAA) objective function. If a target allele
    count (tct) or

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

    # calculate population frequencies (use same data type as tfreq)
    pfreq = sgeno.sum((0,1)) / tfreq.dtype.type(phases)

    # calculate allele availability
    allele_avail = numpy.where(
        tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
        pfreq > 0.0,            # then set TRUE if population has >0 allele freq
        numpy.where(            # else
            tfreq > 0.0,        # if target freq > 0.0
            numpy.logical_and(  # then set TRUE if pop freq is between (0.0,1.0)
                pfreq > 0.0,
                pfreq < 1.0
            ),
            pfreq < 1.0         # else set TRUE if pop freq is less than 1.0
        )
    )

    # calculate population allele availability score using max precision
    paa = (wcoeff * allele_avail).sum(dtype='float64')

    # return score as the specified output data type
    return dtype.type(paa)
