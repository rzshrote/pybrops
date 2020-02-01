import numpy

def stdA_2way(rsel, varA, dtype=numpy.dtype("float64")):
    """
    Given a position vector, return the sum of standard deviations for each
    pairwise cross.

    Bare bones objective function.

    Parameters
    ==========
    rsel : numpy.ndarray
        A 1D array of indices to use for mate pair variance access. The length
        of this array must be divisible by 2!
        Mate are assigned as follows:
            The first element is the male, the second element is the female.
            Example:
                rsel[ 0 ] -> male
                rsel[ 1 ] -> female
                rsel[ 2 ] -> male
                rsel[ 3 ] -> female
                ...
                rsel[n-1] -> male
                rsel[ n ] -> female
    varA : numpy.ndarray
        A 2D array of progeny additive variances. The array should be formatted
        such that the row index is the male index and the col index is the
        female index.
        This function assumes that the varA matrix is completely computed.
    dtype : numpy.dtype
        Function return data type.

    Returns
    =======
    stdsum : numpy.dtype
        A value representing the sum of standard deviations of offspring.
        This is equivalent to the sum of genetic gains assuming constant
        selection intensity and constant heritability.
    """
    # calculate sum of standard deviations from additive var
    stdA = numpy.sqrt(  # take the square root of the additive variance
        varA[           # get a subset of variances
            rsel[0::2], # select 0 to end with step of 2
            rsel[1::2]  # select 1 to end with step of 2
        ]
    ).sum()             # sum up all standard deviations

    # cast as dtype and return stdA
    return dtype.type(stdA)

# lookup table for cross variance objective functions
STDA_DICT = {
    '2-way': stdA_2way
}
