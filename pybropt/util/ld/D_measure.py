import numpy

# TODO: maybe make data return type option, currently only float64
def D_measure(geno):
    """
    Calculate LD as a measure of D as introduced by Lewontin and Kojima (1960).
    This is calculated by t(geno)*geno - t(p)*p. Where geno is a binary matrix
    of allele states and p is a column vector of allele probabilities for the
    dominant state (=1).

    Parameters
    ----------
    geno : numpy.ndarray
        A binary array of allele states. The dtype of 'geno' should be an
        unsigned integer large enough to represent the product M*N.
        Array shape should be (depth, row, column) = (M, N, L) where 'M'
        represents number of chromosome phases, 'N' represents number of
        individuals, 'L' represents number of markers. Array format should be
        the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.

    Returns
    -------
    D : numpy.ndarray
        A 2D matrix of D values for each marker site.
    """
    # get the original shape of the matrix
    s = geno.shape

    # calculate the number of sequences
    nseq = s[0] * s[1]
    nseq_f64 = numpy.float64(nseq)

    # alter the shape of the matrix without copying it.
    # NOTE: assumes C memory layout
    geno.shape = (nseq, s[2])

    # calculate sum of alleles
    allele_count = geno.sum(0)

    # get the allele probabilities
    allele_prob = allele_count / nseq_f64

    # reshape the matrix to be a column matrix.
    allele_prob.shape = (s[2], 1)

    # calculate the matrix
    D = (numpy.matmul(geno.transpose(), geno) / nseq_f64) -
        numpy.matmul(allele_prob, allele_prob.transpose())

    # reshape the geno matrix back to its original shape since it's a pointer.
    geno.shape = s

    return D



def D_measure(geno, rchunk):
    """
    Generator function for calculating an entire LD matrix in chunks.

    Calculate LD as a measure of D as introduced by Lewontin and Kojima (1960).
    This is calculated by t(geno)*geno - t(p)*p. Where geno is a binary matrix
    of allele states and p is a column vector of allele probabilities for the
    dominant state (=1).

    Parameters
    ----------
    geno : numpy.ndarray
        A binary array of allele states. The dtype of 'geno' should be an
        unsigned integer large enough to represent the product M*N.
        Array shape should be (depth, row, column) = (M, N, L) where 'M'
        represents number of chromosome phases, 'N' represents number of
        individuals, 'L' represents number of markers. Array format should be
        the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    rchunk : int
        Number of rows for a calculation chunk. This 

    Returns
    -------
    D : numpy.ndarray
        A 2D matrix of D values for each marker site.
    """
