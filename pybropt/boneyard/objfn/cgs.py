import numpy

def cgs(rslice, geno, coeff):
    """
    Score a population of individuals based on Conventional Genomic Selection
    (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
    Genomic Estimated Breeding Values (GEBV) for a population.

    CGS selects the 'q' individuals with the largest GEBVs.

    Parameters
    ----------
    rslice : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rslice' represents a single individual's row. If 'rslice'
        is None, use all individuals.
        # TODO: performance testing to see which array type is better.
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

    Returns
    -------
    numpy.ndarray
        Returns a 1D array of floating point number representing GEBVs for each
        individual.
    """
    return numpy.dot(       # take dot product of each row in each slice
        geno[:,rslice,:],   # select rows of the matrix
        coeff               # multiply by marker coefficients
    ).sum(0)                # Sum across axis=0 to combine dots products
                            # between slices.
