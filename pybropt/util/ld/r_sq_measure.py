import numpy

def r_sq(rgeno, cgeno, rprob, cprob, r, cycles):
    """
    Calculate r squared linkage disequilibrium matrix.

    Parameters
    ==========
    rsel : numpy.ndarray, tuple, list, None
        A 1D array of indices to use for slicing the 'geno' matrix by row.
        Each index in 'rsel' represents a single individual's row. If 'rsel'
        is None, use all individuals.
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


# def r_sq(rsel, geno, rslice, cslice, r, cycles):
#     """
#     Calculate r squared linkage disequilibrium matrix.
#
#     Parameters
#     ==========
#     rsel : numpy.ndarray, tuple, list, None
#         A 1D array of indices to use for slicing the 'geno' matrix by row.
#         Each index in 'rsel' represents a single individual's row. If 'rsel'
#         is None, use all individuals.
#     geno : numpy.ndarray
#         An array of allele states. The dtype of 'geno' should be 'uint8'. Array
#         shape should be (depth, row, column) = (M, N, L) where 'M' represents
#         number of chromosome phases, 'N' represents number of individuals, 'L'
#         represents number of markers. Array format should be the 'C' format.
#     rslice : slice
#         Row slice. The final r squared matrix has dimensions (L, L). Of this
#         matrix, compute the row slice specified by 'rslice'.
#     cslice : slice
#         Column slice. The final r squared matrix has dimensions (L, L). Of this
#         matrix, compute the column slice specified by 'cslice'.
#     r : numpy.ndarray
#         An array of recombination probabilities. The shape of the array is the
#         same shape specified by 'rslice' and 'cslice': (rslice, cslice).
#     cycles : int, numpy.integer
#         Number of generations of random mating cycles left to the deadline.
#
#     Returns
#     =======
#     r_sq : numpy.ndarray
#         A matrix of size (rslice, cslice) containing r squared values for
#         linkage disequilibrium.
#     """
#     ############################# Compute matrix ##############################
#     # grab a matrix view for target rows and columns
#     row_view = geno[:,rsel,rslice]
#     col_view = geno[:,rsel,cslice]
#
#     # get number of phases
#     phases = row_view.shape[0] * row_view.shape[1]
#
#     # calculate allele proportions for rows and columns
#     # step 1) grab 'rsel' rows; grab columns using 'rslice' or 'cslice'
#     # step 2) sum across depth and rows
#     # step 3) divide by number of phases to get probability of allele
#     row_prob = geno[:,rsel,rslice].sum((0,1)) / phases
#     col_prob = geno[:,rsel,cslice].sum((0,1)) / phases
#
#     # calculate coupling phase matrix
#     coupling = numpy.sum(                       # sum across each phase slice
#         row_view.transpose(0,2,1) @ col_view,   # multiply each phase slice
#         axis=0                                  # phases: axis=0
#     ) / phases                                  # divide by number of phases
#
#     # calculate repulsion phase matrix
#     repulsion = row_prob[:,None] @ col_prob[None,:]
#
#     # calculate D matrix; multiply by (1-r)**cycles for LD decay
#     D = (coupling - repulsion) * ( (1 - r)**cycles )
#
#     # calculate r squared
#     # step 1) square D
#     # step 2) calculate row_prob * (1-row_prob) to get 1D row matrix
#     # step 3) calculate col_prob * (1-col_prob) to get 1D col matrix
#     # step 4) calculate product of row and col probs matrix to get
#     #         (rslice, cslice) dimension matrix.
#     #         This is equivalent to P(A1)*P(1-A1)*P(B1)*P(1-B1)
#     # step 5) divide D_sq by matrix
#     r_sq = (D * D) / (                          # square D
#         (row_prob * (1-row_prob))[:,None] @     # row probs matrix
#         (col_prob * (1-col_prob))[None,:]       # col probs matrix
#     )
#
#     # return r squared metric
#     return r_sq


# TODO: r_sq generator function
def r_sq_gen(rsel, geno, rslice, cslice, r, cycles,
             mem = None, mtype = None, dtype = None):
    return 0

# def r2_measure(phase, cchunk, diag=True):
#     """
#     Generator function for calculating an entire LD matrix in chunks.
#
#     Calculate LD as a measure of D as introduced by Lewontin and Kojima (1960).
#     This is calculated by t(geno)*geno - t(p)*p. Where geno is a binary matrix
#     of allele states and p is a column vector of allele probabilities for the
#     dominant state (=1).
#
#     Parameters
#     ----------
#     phase : numpy.ndarray
#         A binary array of allele states. The dtype of 'phase' should be an
#         unsigned integer large enough to represent the product MN.
#         Array shape should be (row, column) = (MN, L) where 'MN' represents
#         number of chromosome phases (M phases from N individuals), 'L'
#         represents number of markers. Array format should be the 'C' format.
#         NOTE: If you have a 3 dimensional array of size (M, N, L) in a 'C'
#               layout, you can alter phase.shape to press it into a (MN, L)
#               array.
#         TODO: performance testing for 'F' and 'C' array layouts.
#     cchunk : int
#         Number of columns for a calculation chunk. LD matrix chunks will be
#         calculated in chunks of size (cchunk, cchunk) unless the edge of the
#         matrix is encountered, otherwise it will be the remainder size.
#         Suggested chunk size: 512. This corresponds to computational chunks of
#         about 1 Mb with a data size of 4 bytes.
#     diag : boolean
#         If True, only calculate the upper triangle matrix in chunks. This will
#         reduce calculations since an LD matrix is symmetrical. If False,
#         calculate the entire matrix.
#
#     Returns
#     -------
#     D_chunk : numpy.ndarray
#         A 2D matrix chunk of D values for each marker site. Generator order is
#         column major meaning that the calculation frame increases in row index
#         first, then increases in column index. In the case of diag==True, the
#         number of chunks per column set increases until the operation is done.
#     """
#     ##################################################
#     # gather matrix shape data needed for operations #
#     ##################################################
#
#     # get shape as individual
#     rows, columns = geno.shape
#
#     # convert to float64
#     nseq_f64 = numpy.float64(rows)
#
#     ##################################
#     # calculate allele probabilities #
#     ##################################
#
#     # calculate the allele probabilities
#     allele_prob = geno.sum(0) / nseq_f64
#
#     # reshape the matrix to be a column matrix.
#     allele_prob.shape = (columns, 1)
#
#     # calculate A = p(1-p)
#     A = allele_prob * (1 - allele_prob)
#
#     ##############################
#     # calculate LD matrix chunks #
#     ##############################
#
#     ###
#     ### common variabls to diag or non-diag operations
#     ###
#
#     # integer division to get the number of whole chunks
#     num_whole_chunks = columns // cchunk
#
#     # counters for row and column loops, respectively
#     row_pos_cnt = 0
#     col_pos_cnt = 0
#
#     ###
#     ### How this works:
#     ###     Generate chunks of a larger matrix in row major order
#     ### Diagram:
#     ###     +-+-+-+
#     ###     |1|2|3|
#     ###     +-+-+-+
#     ###     |4|5|6|
#     ###     +-+-+-+
#     ###     |7|8|9|
#     ###     +-+-+-+
#     ### Chunks are calulated from 1 to 9 in that order; they are skipped if diag
#     ### variable is True.
#
#     # if we only want the diagonal of the matrix, use a special set of for loops
#     if diag:
#         # for each row chunk in the final matrix
#         for i in range(num_whole_chunks):
#             # select and transpose a section of columns so we get a (c x n)
#             # matrix
#             row_mat = geno[:,row_pos_cnt:(row_pos_cnt + cchunk)].transpose()
#
#             # select a section of the allele probability column vector
#             row_prob = allele_prob[row_pos_cnt:(row_pos_cnt + cchunk),:]
#
#             # select a section of A column vector
#             row_A = A[row_pos_cnt:(row_pos_cnt + cchunk),:]
#
#             # offset the col_pos_cnt variable
#             col_pos_cnt = row_pos_cnt
#
#             # for each column chunk in the final matrix
#             for j in range(i, num_whole_chunks):
#                 # select a section of columns so we get a (n x c) matrix.
#                 col_mat = geno[:,col_pos_cnt:(col_pos_cnt + cchunk)]
#
#                 # select a section of the allele probability column vector
#                 col_prob = allele_prob[col_pos_cnt:(col_pos_cnt + cchunk),:]
#
#                 # select a section of the A column vector
#                 col_A = A[col_pos_cnt:(col_pos_cnt + cchunk),:]
#
#                 # calculate the LD in the chunk
#                 D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                       numpy.matmul(row_prob, col_prob.transpose()))
#
#                 # calculate r_sq
#                 r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#                 # increment col_pos_cnt
#                 col_pos_cnt += cchunk
#
#                 # yield whole thing, otherwise upper triangle matrix.
#                 yield r_sq if row_pos_cnt != col_pos_cnt else numpy.triu(r_sq)
#
#             # calculate a rectangular matrix if one needs to be calculated
#             if columns > col_pos_cnt:
#                 # get everything from col_pos_cnt to end of columns
#                 col_mat = geno[:,col_pos_cnt:columns]
#
#                 # get the last remaining allele probabilities
#                 col_prob = allele_prob[col_pos_cnt:columns,:]
#
#                 # select a section of the A column vector
#                 col_A = A[col_pos_cnt:columns,:]
#
#                 # calculate the LD matrix chunk
#                 D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                       numpy.matmul(row_prob, col_prob.transpose()))
#
#                 # calculate r_sq
#                 r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#                 # yield D. D should never be a square matrix along the diagonal
#                 # it should always be a rectangular matrix.
#                 yield r_sq
#             # increment the counter by cchunk
#             row_pos_cnt += cchunk
#
#         # calculate the last remaining diagonal fi one needs to be calculated
#         if columns > row_pos_cnt:
#             # get everything from row_pos_cnt to the end
#             row_mat = geno[:,row_pos_cnt:columns].transpose()
#
#             # get everything from row_pos_cnt to the end
#             row_prob = allele_prob[row_pos_cnt:columns,:]
#
#             # select a section of the A column vector
#             row_A = A[row_pos_cnt:columns,:]
#
#             # offset the col_pos_cnt variable
#             col_pos_cnt = row_pos_cnt
#
#             # get everything from col_pos_cnt to end of columns
#             col_mat = geno[:,col_pos_cnt:columns]
#
#             # get the last remaining allele probabilities
#             col_prob = allele_prob[col_pos_cnt:columns,:]
#
#             # select a section of the A column vector
#             col_A = A[col_pos_cnt:columns,:]
#
#             # calculate the LD matrix chunk
#             D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                   numpy.matmul(row_prob, col_prob.transpose()))
#
#             # calculate r_sq
#             r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#             # yield upper trangle matrix of D. D should always be a square
#             # matrix along the diagonal.
#             yield numpy.triu(r_sq)
#
#     # else we're doing the entire matrix in chunks
#     else:
#         # for each row chunk in the final matrix
#         for i in range(num_whole_chunks):
#             # select and transpose a section of columns so we get a (c x n)
#             # matrix
#             row_mat = geno[:,row_pos_cnt:(row_pos_cnt + cchunk)].transpose()
#
#             # select a section of the allele probability column vector
#             row_prob = allele_prob[row_pos_cnt:(row_pos_cnt + cchunk),:]
#
#             # select a section of the A column vector
#             row_A = A[row_pos_cnt:(row_pos_cnt + cchunk),:]
#
#             # offset the col_pos_cnt variable
#             col_pos_cnt = 0
#
#             # for each column chunk in the final matrix
#             for j in range(num_whole_chunks):
#                 # select a section of columns so we get a (n x c) matrix.
#                 col_mat = geno[:,col_pos_cnt:(col_pos_cnt + cchunk)]
#
#                 # select a section of the allele probability column vector
#                 col_prob = allele_prob[col_pos_cnt:(col_pos_cnt + cchunk),:]
#
#                 # select a section of the A column vector
#                 col_A = A[col_pos_cnt:(col_pos_cnt + cchunk),:]
#
#                 # calculate the LD in the chunk
#                 D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                       numpy.matmul(row_prob, col_prob.transpose()))
#
#                 # increment col_pos_cnt
#                 col_pos_cnt += cchunk
#
#                 # calculate r_sq
#                 r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#                 # yield r squared chunk
#                 yield r_sq
#
#             # calculate a rectangular matrix if one needs to be calculated
#             if columns > col_pos_cnt:
#                 # get everything from col_pos_cnt to end of columns
#                 col_mat = geno[:,col_pos_cnt:columns]
#
#                 # get the last remaining allele probabilities
#                 col_prob = allele_prob[col_pos_cnt:columns,:]
#
#                 # select a section of the A column vector
#                 col_A = A[col_pos_cnt:columns,:]
#
#                 # calculate the LD matrix chunk
#                 D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                       numpy.matmul(row_prob, col_prob.transpose()))
#
#                 # calculate r_sq
#                 r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#                 # yield r_sq. r_sq should never be a square matrix along the
#                 # diagonal it should always be a rectangular matrix.
#                 yield D
#             # increment the counter by cchunk
#             row_pos_cnt += cchunk
#
#         # calculate the last remaining diagonal fi one needs to be calculated
#         if columns > row_pos_cnt:
#             # get everything from row_pos_cnt to the end
#             row_mat = geno[:,row_pos_cnt:columns].transpose()
#
#             # get everything from row_pos_cnt to the end
#             row_prob = allele_prob[row_pos_cnt:columns,:]
#
#             # select a section of the A column vector
#             row_A = A[row_pos_cnt:columns,:]
#
#             # offset the col_pos_cnt variable
#             col_pos_cnt = 0
#
#             # for each column chunk in the final matrix
#             for j in range(num_whole_chunks):
#                 # select a section of columns so we get a (n x c) matrix.
#                 col_mat = geno[:,col_pos_cnt:(col_pos_cnt + cchunk)]
#
#                 # select a section of the allele probability column vector
#                 col_prob = allele_prob[col_pos_cnt:(col_pos_cnt + cchunk),:]
#
#                 # select a section of the A column vector
#                 col_A = A[col_pos_cnt:(col_pos_cnt + cchunk),:]
#
#                 # calculate the LD in the chunk
#                 D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                       numpy.matmul(row_prob, col_prob.transpose()))
#
#                 # calculate r_sq
#                 r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#                 # increment col_pos_cnt
#                 col_pos_cnt += cchunk
#
#                 # yield whole thing, otherwise upper triangle matrix.
#                 yield D
#
#             # calculate a rectangular matrix if one needs to be calculated
#             if columns > col_pos_cnt:
#                 # get everything from col_pos_cnt to end of columns
#                 col_mat = geno[:,col_pos_cnt:columns]
#
#                 # get the last remaining allele probabilities
#                 col_prob = allele_prob[col_pos_cnt:columns,:]
#
#                 # select a section of the A column vector
#                 col_A = A[col_pos_cnt:columns,:]
#
#                 # calculate the LD matrix chunk
#                 D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
#                       numpy.matmul(row_prob, col_prob.transpose()))
#
#                 # calculate r_sq
#                 r_sq = (D*D) / numpy.matmul(row_A, col_A.transpose())
#
#                 # yield r squared
#                 yield r_sq
