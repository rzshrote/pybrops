import numpy

# TODO: clean me!

# good function: can be used with optimzation algorithms
def opv(rslice, hcoeff):
    """
    Score a population of individuals based on Optimal Population Value (OPV)
    (Goiffon et al., 2017). Scoring for OPV is defined as the Genomic Estimated
    Breeding Value (GEBV) of the best inbred progeny produced from a breeding
    population allowing for an infinite number of intercrossings and
    generations. Haplotypes can be specified in OPV and markers within the
    haplotypes are assumed to always segregate with the haplotype. In simpler
    terms, OPV collects all of the best haplotype blocks and calculates the GEBV
    for those haplotypes for a given set of individuals.

    OPV selects the 'q' individuals to maximize the maximum GEBV possible within
    a population consisting of the 'q' selected individuals.

    Parameters
    ----------
    hcoeff : numpy.ndarray
        An array of coefficients for haplotype effects. The dtype of 'hcoeff'
        should be either 'float32' or 'float64'. Array shape should be
        (depth, row, column) = (M, N, H) where 'M' represents number of
        chromosome phases, 'N' represents number of individuals, 'H' represents
        number of haplotypes. Array format should be the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    rslice : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rslice' represents a single individual's row. If 'rslice'
        is None, use all individuals.
        # TODO: performance testing to see which array type is better.

    Returns
    -------
    float
        Returns a floating point number representing the OPV score. Data type
        depends on the dtype of 'coeff'.
    """
    ############################################################################
    # Algorithm:
    #     Take the array max (numpy.amax) along the chromosome and individual
    #     axes (axis=(0,1), respectively), then take the sum of the 1D array to
    #     arrive at OPV GEBV.
    ############################################################################

    return numpy.amax(      # array max function
        hcoeff[:,rslice,:], # subset hcoeff rows
        axis=(0,1)          # calculate max along slice and row axes
                            # NOTE: this is on a per-haplotype basis.
    ).sum()                 # take sum of maximum haplotypes


################################################################################
# Calculate OPV using a haplotype value matrix.
################################################################################
def opv_hcoeff(hcoeff, rslice = None):
    """
    Score a population of individuals based on Optimal Population Value (OPV)
    (Goiffon et al., 2017). Scoring for OPV is defined as the Genomic Estimated
    Breeding Value (GEBV) of the best inbred progeny produced from a breeding
    population allowing for an infinite number of intercrossings and
    generations. Haplotypes can be specified in OPV and markers within the
    haplotypes are assumed to always segregate with the haplotype. In simpler
    terms, OPV collects all of the best haplotype blocks and calculates the GEBV
    for those haplotypes for a given set of individuals.

    OPV selects the 'q' individuals to maximize the maximum GEBV possible within
    a population consisting of the 'q' selected individuals.

    Parameters
    ----------
    hcoeff : numpy.ndarray
        An array of coefficients for haplotype effects. The dtype of 'hcoeff'
        should be either 'float32' or 'float64'. Array shape should be
        (depth, row, column) = (M, N, H) where 'M' represents number of
        chromosome phases, 'N' represents number of individuals, 'H' represents
        number of haplotypes. Array format should be the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    rslice : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rslice' represents a single individual's row. If 'rslice'
        is None, use all individuals.
        # TODO: performance testing to see which array type is better.

    Returns
    -------
    float
        Returns a floating point number representing the OPV score. Data type
        depends on the dtype of 'coeff'.
    """
    ############################################################################
    # Algorithm:
    #     Take the array max (numpy.amax) along the chromosome and individual
    #     axes (axis=(0,1), respectively), then take the sum of the 1D array to
    #     arrive at OPV GEBV.
    ############################################################################

    # if we only want specific rows, subset those rows, then calculate sum.
    if rslice != None:
        return numpy.amax(      # array max function
            hcoeff[:,rslice,:], # subset hcoeff rows
            axis=(0,1)          # calculate max along slice and row axes
                                # NOTE: this is on a per-haplotype basis.
        ).sum()                 # take sum of maximum haplotypes

    # otherwise, we don't care and consider all individuals
    return numpy.amax(          # array max function
        hcoeff,                 # all hcoeff rows
        axis=(0,1)              # calculate max along slice and row axes
                                # NOTE: this is on a per-haplotype basis.
    ).sum()                     # take sum of maximum haplotypes

################################################################################
# Calculate OPV using genotype, coefficient, hbin matrices.
# I will probably depricate this.
################################################################################
def opv_geno(geno, coeff, hbin = None, rslice = None):
    """
    Score a population of individuals based on Optimal Population Value (OPV)
    (Goiffon et al., 2017). Scoring for OPV is defined as the Genomic Estimated
    Breeding Value (GEBV) of the best inbred progeny produced from a breeding
    population allowing for an infinite number of intercrossings and
    generations. Haplotypes can be specified in OPV and markers within the
    haplotypes are assumed to always segregate with the haplotype. In simpler
    terms, OPV collects all of the best haplotype blocks and calculates the GEBV
    for those haplotypes for a given set of individuals.

    OPV selects the 'q' individuals to maximize the maximum GEBV possible within
    a population consisting of the 'q' selected individuals.

    Parameters
    ----------
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    coeff : numpy.ndarray
        An array of coefficients for allele effects. The dtype of 'coeff'
        should be either 'float32' or 'float64'. This array should be single
        dimensional (column,) = (N,) where 'N' represents number of markers.
    hbin : numpy.ndarray
        An array of haplotype groupings. This array should be single dimensional
        (column,) = (N,) where 'N' represents number of markers. Groupings
        should be specified by different integers.
        Example:
            hbin = array([0, 0, 1, 1, 2, 2], dtype='uint32')
            The array within 'hbin' will group indices (0,1), (2,3), (4,5).
            This will group markers 1-2, 3-4, 5-6 into three respective
            haplotypes.
    rslice : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rslice' represents a single individual's row. If 'rslice'
        is None, use all individuals.
        # TODO: performance testing to see which array type is better.

    Returns
    -------
    float
        Returns a floating point number representing the OPV score. Data type
        depends on the dtype of 'coeff'.
    """
    ############################################################################
    # Algorithm:
    #     Multiply the 'geno' matrix by the 'coeff' matrix to get values of each
    #     marker, then bin markers into their respective haplotypes.
    ############################################################################

    # declare haplotype coefficient matrix variable (null)
    hcoeff = None

    # if hbin is not specified, then each marker is a haplotype.
    if hbin == None:

        # if specific rows are not specified, then use all of them.
        if rslice == None:
            # WARNING: Current implementation assumes 'C' matrix layout.
            # Mulitply genotype by coefficients; store into hcoeff
            hcoeff = numpy.multiply(geno, coeff, order='C')

        # else specific rows have been specified
        else:
            # WARNING: Current implementation assumes 'C' matrix layout.
            # Multply select rows by coefficients; store into hcoeff
            hcoeff = numpy.multiply(geno[:,rslice,:], coeff, order='C')

    # else hbin is specified and we need to group into specific haplotype bins.
    else:
        # get the dimensions of the array
        slices, rows, columns = ar3d_C.shape

        # make a null variable for riter; this will be used successively in for
        # loops below.
        riter = None

        # allocate memory for empty array using slice, row, haplotype counts for
        # dimensions and the coeff dtype as the dtype.
        hcoeff = np.empty((slices, rows, np.amax()+1), dtype=coeff.dtype)

        # for each slice
        for slice in range(slices):
            # make a row iterator from range() if no rslice rows are specified,
            # else make iterator from the rslice array.
            riter = range(rows) if rslice == None else numpy.nditer(rslice)

            # for each selected row
            for row in riter:
                # copy weighted bincounts to the array
                hcoeff[slice,row,:] = numpy.bincount( # bincount function
                    hbin,                             # bin group specifications
                    weights = numpy.multiply(         # bin weights from product
                        geno[slice,row,:],            # row genotype data
                        coeff                         # coeff's for each column
                    )
                )



    ############################################################################
    # Algorithm:
    #     Take the array max (numpy.amax) along the chromosome and individual
    #     axes (axis=(0,1), respectively), then take the sum of the 1D array to
    #     arrive at OPV GEBV.
    ############################################################################

    # if we only want specific rows, subset those rows
    if rslice != None:
        return numpy.amax(hcoeff[:,rslice,:], axis=(0,1)).sum()

    # otherwise, we don't care and consider all individuals
    return numpy.amax(hcoeff, axis=(0,1)).sum()
