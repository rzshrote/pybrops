import numpy

def pau(rsel, geno, wcoeff, tfreq, dtype=numpy.dtype("float64")):
    """
    Population Allele Unvailability (PAU) objective function.
        The goal is to minimize this function. Lower is better.
        This is a bare bones function. Minimal error checking.

    Given the provided genotype matrix 'geno' and row selections from it 'rsel',
    calculate the selection allele freq. From the selection allele frequencies
    and the target allele frequencies, determine if the target frequencies
    cannot be attained after unlimited generations and selection rounds.
    Multiply this vector by a weight coefficients vector 'wcoeff'.

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
    dtype : numpy.dtype
        Specify the return type for this function.

    Returns
    =======
    pau : numpy.dtype
        The PAU score of the dtype specified by 'dtype'.
    """
    # if 'wcoeff' matrix does not have a dtype of float64, convert it.
    # use numpy.float64(a) instead of a.astype('float64') for speed
    if wcoeff.dtype != 'float64':
        wcoeff = numpy.float64(wcoeff)

    # generate a view of the geno matrix that only contains 'rsel' rows.
    sgeno = geno[:,rsel,:]

    # calculate the number of chromosome phases present (depth * rows)
    phases = sgeno.shape[0] * sgeno.shape[1]

    # calculate population frequencies (use float64)
    pfreq = sgeno.sum((0,1)) / numpy.float64(phases)

    # calculate some inequalities for use multiple times
    pfreq_lteq_0 = (pfreq <= 0.0)   # is population freq <= 0.0
    pfreq_gteq_1 = (pfreq >= 1.0)   # is population freq >= 1.0

    # calculate allele unavailability
    allele_unavail = numpy.where(
        tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
        pfreq_lteq_0,           # then set TRUE if population has =0 allele freq
        numpy.where(            # else
            tfreq > 0.0,        # if 0.0 < target freq < 1.0
            numpy.logical_or(   # then set TRUE if pop freq is outside (0.0,1.0)
                pfreq_lteq_0,
                pfreq_gteq_1
            ),
            pfreq_gteq_1        # else set TRUE if pop freq is >= 1.0
        )
    )

    # calculate population allele unavailability score; will upgrade to float64
    pau = numpy.dot(wcoeff, allele_unavail)

    # return score as the specified output data type
    return dtype.type(pau)
