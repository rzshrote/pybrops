def hcoeff(geno, coeff, hbin = None, rslice = None):
    """
    Group markers together into haplotype blocks and compute a haplotype matrix.
    Each value within the matrix represents the coefficient of a haplotype.

    Parameters
    ----------
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosomes, 'N' represents number of individuals, 'L'
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
    numpy.ndarray
        Returns an array of floating point numbers representing the coefficients
        of haplotypes. Array data type depends on the dtype of 'coeff'.
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

    # return our haplotypes
    return hcoeff
