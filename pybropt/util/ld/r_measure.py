import numpy

def r_(rsel, geno, rslice, cslice, r, cycles, dtype = None):
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
    r_ : numpy.ndarray
        A matrix of size (rslice, cslice) containing r values for linkage
        disequilibrium.
    """
    ############################## Parse 'dtype' ##############################
    if dtype is None:
        dtype = r.dtype

    ############################# Compute matrix ##############################
    # get number of phases
    phases = geno.shape[0] * geno.shape[1]
    depth, rows, columns = geno.shape

    # grab a matrix view for target rows and columns
    row_view = geno[:,rsel,rslice]
    col_view = geno[:,rsel,cslice]

    # calculate allele proportions for rows and columns
    # step 1) grab 'rsel' rows; grab columns using 'rslice' or 'cslice'
    # step 2) sum across depth and rows
    # step 3) divide by number of phases to get probability of allele
    row_prob = geno[:,rsel,rslice].sum((0,1)) / phases
    col_prob = geno[:,rsel,cslice].sum((0,1)) / phases

    # calculate coupling phase matrix
    coupling = numpy.sum(                       # sum across each phase slice
        row_view.transpose(0,2,1) @ col_view,   # multiply each phase slice
        axis=0                                  # phases: axis=0
    ) / phases                                  # divide by number of phases

    # calculate repulsion phase matrix
    repulsion = row_prob[:,None] @ col_prob[None,:]

    # calculate D matrix; multiply by (1-r)**cycles for LD decay
    D = (coupling - repulsion) * ( (1 - r)**cycles )

    # calculate r
    r_ = D / numpy.sqrt(                        # calculate sqrt(r_sq)
        (row_prob * (1-row_prob))[:,None] @     # row probs matrix
        (col_prob * (1-col_prob))[None,:]       # col probs matrix
    )

    # return r metric
    return r_


# TODO: r_sq generator function
def r_sq_gen(rsel, geno, rslice, cslice, r, cycles,
             mem = None, mtype = None, dtype = None):
    return 0
