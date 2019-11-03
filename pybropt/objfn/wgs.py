import numpy

def wgs(geno, coeff, weight, rslice = None):
    """
    Score a population of individuals based on Weighted Genomic Selection (WGS)
    (Goddard, 2009; Jannink, 2010). Scoring for WGS is defined as the adjusted
    Genomic Estimated Breeding Values (GEBV) for an individual. Marker
    coefficients are adjusted according to the frequency of the most beneficial
    allele.

    WGS selects the 'q' individuals with the largest weighted GEBVs.

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
    weight : numpy.ndarray
        An array of weights for allele effect coefficients. The dtype of
        'weight' should be either 'float32' or 'float64'. This array should be
        single dimensional (column,) = (N,) where 'N' represents number of
        markers.
    rslice : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rslice' represents a single individual's row. If 'rslice'
        is None, use all individuals.
        # TODO: performance testing to see which array type is better.

    Returns
    -------
    numpy.ndarray
        Returns a 1D array of floating point number representing GEBVs for each
        individual.
    """
    if rslice != None:
        return numpy.dot(geno[:,rslice,:], coeff / numpy.sqrt(weight)).sum(0)
    return numpy.dot(geno, coeff / numpy.sqrt(weight)).sum(0)
