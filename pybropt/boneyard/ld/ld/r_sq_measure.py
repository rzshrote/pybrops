import numpy

def r_sq(rgeno, cgeno, rprob, cprob, r, cycles):
    """
    Calculate r squared linkage disequilibrium matrix.

    Parameters
    ==========
    rgeno : numpy.ndarray
        A binary array of allele states for output matrix rows. This array
        determines row number of the output matrix. The dtype of 'rgeno'
        should be 'uint8'. Array shape should be (depth, row, column) =
        (M, N, L) where 'M' represents number of chromosome phases, 'N'
        represents number of individuals, 'L' represents number of markers.
        Dimensions of 'rgeno' and 'cgeno' should be compatible. This means:
            rgeno.M == cgeno.M
            rgeno.N == cgeno.N
        rgeno.L and cgeno.L do not have to be the same size:
            rgeno.L determines output matrix row number
            cgeno.L determines output matrix col number
        Array format should be the 'C' format.
    cgeno : numpy.ndarray
        A binary array of allele states for output matrix columns. This array
        determines column number of the output matrix. The dtype of 'cgeno'
        should be 'uint8'. Array shape should be (depth, row, column) =
        (M, N, L) where 'M' represents number of chromosome phases, 'N'
        represents number of individuals, 'L' represents number of markers.
        Dimensions of 'rgeno' and 'cgeno' should be compatible. This means:
            rgeno.M == cgeno.M
            rgeno.N == cgeno.N
        rgeno.L and cgeno.L do not have to be the same size:
            rgeno.L determines output matrix row number
            cgeno.L determines output matrix col number
        Array format should be the 'C' format.
    rprob : numpy.ndarray
        An array of allele probabilities for output matrix rows. This array
        determines row number of the output matrix. The length of 'rprob'
        should be equivalent to the column number in 'rgeno'.
    cprob : numpy.ndarray
        An array of allele probabilities for output matrix columns. This array
        determines column number of the output matrix. The length of 'cprob'
        should be equivalent to the column number in 'cgeno'.
    r : numpy.ndarray
        An array of recombination probabilities. The shape of the array is the
        same shape specified by the lengths of 'rprob' and 'cprob'.
        r.shape = ( rows=len(rprob), cols=len(cprob) )
    cycles : int, numpy.integer
        Number of generations of random mating cycles to apply after computing
        the linkage disequilibrium.

    Returns
    =======
    r_sq : numpy.ndarray
        A matrix of size (rslice, cslice) containing r squared values for
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

    # calculate r squared
    # step 1) square D
    # step 2) calculate rprob * (1-rprob) to get 1D row matrix
    # step 3) calculate cprob * (1-cprob) to get 1D col matrix
    # step 4) calculate product of 'rprob' and 'cprob' to get
    #         (len(rprob), len(cprob)) dimension matrix.
    #         This is equivalent to P(A1)*P(1-A1)*P(B1)*P(1-B1)
    # step 5) divide D_sq by matrix
    r_sq = (D * D) / (                    # square D
        (rprob * (1-rprob))[:,None] @     # row probs matrix
        (cprob * (1-cprob))[None,:]       # col probs matrix
    )

    # return r squared metric
    return r_sq


def r_sq_imp(rgeno, cgeno, rprob, cprob, ravail, cavail, r, generations):
    """
    Calculate r squared linkage disequilibrium matrix, but imputing NaN's with
    pairwise availablility matrices.

    Parameters
    ==========
    rgeno : numpy.ndarray
        A binary array of allele states for output matrix rows. This array
        determines row number of the output matrix. The dtype of 'rgeno'
        should be 'uint8'. Array shape should be (depth, row, column) =
        (M, N, L) where 'M' represents number of chromosome phases, 'N'
        represents number of individuals, 'L' represents number of markers.
        Dimensions of 'rgeno' and 'cgeno' should be compatible. This means:
            rgeno.M == cgeno.M
            rgeno.N == cgeno.N
        rgeno.L and cgeno.L do not have to be the same size:
            rgeno.L determines output matrix row number
            cgeno.L determines output matrix col number
        Array format should be the 'C' format.
    cgeno : numpy.ndarray
        A binary array of allele states for output matrix columns. This array
        determines column number of the output matrix. The dtype of 'cgeno'
        should be 'uint8'. Array shape should be (depth, row, column) =
        (M, N, L) where 'M' represents number of chromosome phases, 'N'
        represents number of individuals, 'L' represents number of markers.
        Dimensions of 'rgeno' and 'cgeno' should be compatible. This means:
            rgeno.M == cgeno.M
            rgeno.N == cgeno.N
        rgeno.L and cgeno.L do not have to be the same size:
            rgeno.L determines output matrix row number
            cgeno.L determines output matrix col number
        Array format should be the 'C' format.
    rprob : numpy.ndarray
        An array of allele probabilities for output matrix rows. This array
        determines row number of the output matrix. The length of 'rprob'
        should be equivalent to the column number in 'rgeno'.
    cprob : numpy.ndarray
        An array of allele probabilities for output matrix columns. This array
        determines column number of the output matrix. The length of 'cprob'
        should be equivalent to the column number in 'cgeno'.
    ravail : numpy.ndarray
        An array of numpy.bool's that indicate if the desired allele(s) is(are)
        available at the row site of interest. This is used to impute LD
        metrics in the event that division by zero occurs (NaN).
    cavail : numpy.ndarray
        An array of numpy.bool's that indicate if the desired allele(s) is(are)
        available at the row site of interest. This is used to impute LD
        metrics in the event that division by zero occurs (NaN).
    r : numpy.ndarray
        An array of recombination probabilities. The shape of the array is the
        same shape specified by the lengths of 'rprob' and 'cprob'.
        r.shape = ( rows=len(rprob), cols=len(cprob) )
    generations : int, numpy.integer
        Number of generations of random mating cycles to apply after computing
        the linkage disequilibrium.

    Returns
    =======
    r_sq : numpy.ndarray
        A matrix of size (rslice, cslice) containing r squared values for
        linkage disequilibrium.
    """
    # convert 'rprob' and 'cprob' to float64 if needed
    if rprob.dtype != 'float64':        # if rprob is not float64
        rprob = numpy.float64(rprob)    # cast to float64
    if cprob.dtype != 'float64':        # if cprob is not float64
        cprob = numpy.float64(cprob)    # cast to float64

    ############################# Compute matrix ##############################
    # get number of phases from row matrix
    phases = rgeno.shape[0] * rgeno.shape[1]

    # calculate coupling phase matrix
    coupling = numpy.sum(                 # sum across each phase slice
        rgeno.transpose(0,2,1) @ cgeno,   # multiply each phase slice
        axis=0                            # sum across phases: axis=0
    ) / numpy.float64(phases)             # divide by phases; use float64

    # calculate repulsion phase matrix
    repulsion = numpy.outer(rprob, cprob)

    ############# r_sq matrix assembly #############
    # NOTE: r_sq is modified in place to reduce memory space
    #       It contains different values at different time points
    #       Step 1) r_sq contains D matrix
    #       Step 2) r_sq contains D_sq matrix
    #       Step 4) r_sq contains r_sq and NaN
    #       Step 7) r_sq contains r_sq and imputation

    # Step 1) calculate D matrix; multiply by (1-r)**cycles for LD decay
    # since coupling is guaranteed to be float64, this is a float64 matrix
    r_sq = (coupling - repulsion) * ( (1 - r)**generations )

    # Step 2) now square D to get D_sq
    r_sq *= r_sq

    # Step 3) calculate P(A1)*P(1-A1)*P(B1)*P(1-B1) outer product
    denom = numpy.outer(            # this will be denominator in r_sq equation
        (rprob * (1.0 - rprob)),    # row probs matrix
        (cprob * (1.0 - cprob))     # col probs matrix
    )

    # Step 4) divide by D_sq by denom, this will generate NaN's and generate
    # a warning. We do this to take advantage of vectorized math.
    r_sq /= denom

    # Step 5) calculate availablility outer product for imputation purposes
    avail = numpy.outer(ravail, cavail)

    # Step 6) make 0.0 mask for denom
    mask_zero = (denom == 0.0)

    # Step 7) for zero denominators, set to pairwise allele availablility
    r_sq[mask_zero] = avail[mask_zero]

    # return r squared metric
    return r_sq



# TODO: r_sq generator function
def r_sq_gen(rsel, geno, rslice, cslice, r, cycles,
             mem = None, mtype = None, dtype = None):
    return 0
