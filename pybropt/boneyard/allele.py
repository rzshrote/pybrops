import numpy

def allele_freq(rslice, geno):
    """
    Calculate allele frequencies across the genome

    Parameters
    ==========
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

    Returns
    =======
    numpy.ndarray
        Returns a 1D array of floating point number representing allele
        frequencies for each locus.
    """

    return geno[:,rslice,:].sum(0,1) / geno[:,rslice,:].shape[1]
