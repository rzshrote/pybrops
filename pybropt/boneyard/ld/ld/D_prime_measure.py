import numpy

def D_prime(rgeno, cgeno, rprob, cprob, r, cycles):
    """
    Calculate r squared linkage disequilibrium matrix.

    Parameters
    ==========
    rsel : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing the 'geno' matrix by row.
        Each index in 'rsel' represents a single individual's row. If 'rsel'
        is None, use all individuals.
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
    rslice : slice
        Row slice. The final r squared matrix has dimensions (L, L). Of this
        matrix, compute the row slice specified by 'rslice'.
    cslice : slice
        Column slice. The final r squared matrix has dimensions (L, L). Of this
        matrix, compute the column slice specified by 'cslice'.
    r : numpy.ndarray
        An array of recombination probabilities. The shape of the array is the
        same shape specified by 'rslice' and 'cslice': (rslice, cslice).
    cycles : int, numpy.integer
        Number of generations of random mating cycles left to the deadline.
    dtype : numpy.dtype
        The type of the output array. If dtype is not given, infer the data
        type from the 'r' argument.

    Returns
    =======
    D_prime : numpy.ndarray
        A matrix of size (rslice, cslice) containing D prime values for
        linkage disequilibrium.
    """
    ############################# Compute matrix ##############################
    # get number of phases from row matrix
    phases = rgeno.shape[0] * rgeno.shape[1]

    # calculate coupling phase matrix
    coupling = numpy.sum(                 # sum across each phase slice
        rgeno.transpose(0,2,1) @ cgeno,   # multiply each phase slice
        axis=0                            # phases: axis=0
    ) / rprob.dtype.type(phases)          # divide by phases; use rprob dtype

    # calculate repulsion phase matrix
    repulsion = rprob[:,None] @ cprob[None,:]

    # calculate D matrix; multiply by (1-r)**cycles for LD decay
    D = (coupling - repulsion) * ( (1 - r)**cycles )

    # calculate D max
    D_max = numpy.where(
        D < 0,          # if D < 0
        numpy.maximum(  # then take max( -P(A1)P(B1), -(1-P(A1))(1-P(B1)) )
            -rprob[:,None] @ cprob[None,:],
            -(1-rprob)[:,None] @ (1-cprob)[None,:]
        ),
        numpy.minimum(  # else take min( P(A1)(1-P(B1)), (1-P(A1))P(B1) )
            rprob[:,None] @ (1-cprob)[None,:],
            (1-rprob)[:,None] @ cprob[None,:]
        )
    )

    # calculate D prime
    D_prime = D / D_max

    # return D prime
    return D_prime



# TODO: D_prime generator function
def D_prime_gen(rsel, geno, rslice, cslice, r, cycles,
             mem = None, mtype = None, dtype = None):
    return 0
