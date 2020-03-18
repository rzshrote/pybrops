import numpy

# TODO: maybe make data return type option, currently only float64
def D_measure(phase):
    """
    Calculate LD as a measure of D as introduced by Lewontin and Kojima (1960).
    This is calculated by t(geno)*geno - t(p)*p. Where geno is a binary matrix
    of allele states and p is a column vector of allele probabilities for the
    dominant state (=1).

    Parameters
    ----------
    phase : numpy.ndarray
        A binary array of allele states. The dtype of 'phase' should be an
        unsigned integer large enough to represent the product MN.
        Array shape should be (row, column) = (MN, L) where 'MN' represents
        number of chromosome phases (M phases from N individuals), 'L'
        represents number of markers. Array format should be the 'C' format.
        NOTE: If you have a 3 dimensional array of size (M, N, L) in a 'C'
              layout, you can alter phase.shape to press it into a (MN, L)
              array.
        TODO: performance testing for 'F' and 'C' array layouts.

    Returns
    -------
    D : numpy.ndarray
        A 2D matrix of D values for each marker site.
    """
    # get the dimensions of the matrix
    rows, columns = geno.shape

    # convert the number of sequences (rows) to float64
    nseq_f64 = numpy.float64(rows)

    # calculate allele sums (axis=0), divide by nseq_f64 to get probabilities
    allele_prob = geno.sum(0) / nseq_f64

    # reshape the matrix to be a column matrix.
    allele_prob.shape = (columns, 1)

    # calculate the matrix and return it
    return ((numpy.matmul(geno.transpose(), geno) / nseq_f64) -
             numpy.matmul(allele_prob, allele_prob.transpose()))



def D_measure(phase, cchunk, diag=True):
    """
    Generator function for calculating an entire LD matrix in chunks.

    Calculate LD as a measure of D as introduced by Lewontin and Kojima (1960).
    This is calculated by t(geno)*geno - t(p)*p. Where geno is a binary matrix
    of allele states and p is a column vector of allele probabilities for the
    dominant state (=1).

    Parameters
    ----------
    phase : numpy.ndarray
        A binary array of allele states. The dtype of 'phase' should be an
        unsigned integer large enough to represent the product MN.
        Array shape should be (row, column) = (MN, L) where 'MN' represents
        number of chromosome phases (M phases from N individuals), 'L'
        represents number of markers. Array format should be the 'C' format.
        NOTE: If you have a 3 dimensional array of size (M, N, L) in a 'C'
              layout, you can alter phase.shape to press it into a (MN, L)
              array.
        TODO: performance testing for 'F' and 'C' array layouts.
    cchunk : int
        Number of columns for a calculation chunk. LD matrix chunks will be
        calculated in chunks of size (cchunk, cchunk) unless the edge of the
        matrix is encountered, otherwise it will be the remainder size.
        Suggested chunk size: 512. This corresponds to computational chunks of
        about 1 Mb with a data size of 4 bytes.
    diag : boolean
        If True, only calculate the upper triangle matrix in chunks. This will
        reduce calculations since an LD matrix is symmetrical. If False,
        calculate the entire matrix.

    Returns
    -------
    D_chunk : numpy.ndarray
        A 2D matrix chunk of D values for each marker site. Generator order is
        column major meaning that the calculation frame increases in row index
        first, then increases in column index. In the case of diag==True, the
        number of chunks per column set increases until the operation is done.
    """
    ##################################################
    # gather matrix shape data needed for operations #
    ##################################################

    # get shape as individual
    rows, columns = geno.shape

    # convert to float64
    nseq_f64 = numpy.float64(rows)

    ##################################
    # calculate allele probabilities #
    ##################################

    # calculate the allele probabilities
    allele_prob = geno.sum(0) / nseq_f64

    # reshape the matrix to be a column matrix.
    allele_prob.shape = (columns, 1)

    ##############################
    # calculate LD matrix chunks #
    ##############################

    ###
    ### common variabls to diag or non-diag operations
    ###

    # integer division to get the number of whole chunks
    num_whole_chunks = columns // cchunk

    # counters for row and column loops, respectively
    row_pos_cnt = 0
    col_pos_cnt = 0

    ###
    ### How this works:
    ###     Generate chunks of a larger matrix in row major order
    ### Diagram:
    ###     +-+-+-+
    ###     |1|2|3|
    ###     +-+-+-+
    ###     |4|5|6|
    ###     +-+-+-+
    ###     |7|8|9|
    ###     +-+-+-+
    ### Chunks are calulated from 1 to 9 in that order; they are skipped if diag
    ### variable is True.

    # if we only want the diagonal of the matrix, use a special set of for loops
    if diag:
        # for each row chunk in the final matrix
        for i in range(num_whole_chunks):
            # select and transpose a section of columns so we get a (c x n)
            # matrix
            row_mat = geno[:,row_pos_cnt:(row_pos_cnt + cchunk)].transpose()

            # select a section of the allele probability column vector
            row_prob = allele_prob[row_pos_cnt:(row_pos_cnt + cchunk),:]

            # offset the col_pos_cnt variable
            col_pos_cnt = row_pos_cnt

            # for each column chunk in the final matrix
            for j in range(i, num_whole_chunks):
                # select a section of columns so we get a (n x c) matrix.
                col_mat = geno[:,col_pos_cnt:(col_pos_cnt + cchunk)]

                # select a section of the allele probability column vector
                col_prob = allele_prob[col_pos_cnt:(col_pos_cnt + cchunk),:]

                # calculate the LD in the chunk
                D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                      numpy.matmul(row_prob, col_prob.transpose()))

                # increment col_pos_cnt
                col_pos_cnt += cchunk

                # yield whole thing, otherwise upper triangle matrix.
                yield D if row_pos_cnt != col_pos_cnt else numpy.triu(D)

            # calculate a rectangular matrix if one needs to be calculated
            if columns > col_pos_cnt:
                # get everything from col_pos_cnt to end of columns
                col_mat = geno[:,col_pos_cnt:columns]

                # get the last remaining allele probabilities
                col_prob = allele_prob[col_pos_cnt:columns,:]

                # calculate the LD matrix chunk
                D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                      numpy.matmul(row_prob, col_prob.transpose()))

                # yield D. D should never be a square matrix along the diagonal
                # it should always be a rectangular matrix.
                yield D
            # increment the counter by cchunk
            row_pos_cnt += cchunk

        # calculate the last remaining diagonal fi one needs to be calculated
        if columns > row_pos_cnt:
            # get everything from row_pos_cnt to the end
            row_mat = geno[:,row_pos_cnt:columns].transpose()

            # get everything from row_pos_cnt to the end
            row_prob = allele_prob[row_pos_cnt:columns,:]

            # offset the col_pos_cnt variable
            col_pos_cnt = row_pos_cnt

            # get everything from col_pos_cnt to end of columns
            col_mat = geno[:,col_pos_cnt:columns]

            # get the last remaining allele probabilities
            col_prob = allele_prob[col_pos_cnt:columns,:]

            # calculate the LD matrix chunk
            D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                  numpy.matmul(row_prob, col_prob.transpose()))

            # yield upper trangle matrix of D. D should always be a square
            # matrix along the diagonal.
            yield numpy.triu(D)

    # else we're doing the entire matrix in chunks
    else:
        # for each row chunk in the final matrix
        for i in range(num_whole_chunks):
            # select and transpose a section of columns so we get a (c x n)
            # matrix
            row_mat = geno[:,row_pos_cnt:(row_pos_cnt + cchunk)].transpose()

            # select a section of the allele probability column vector
            row_prob = allele_prob[row_pos_cnt:(row_pos_cnt + cchunk),:]

            # offset the col_pos_cnt variable
            col_pos_cnt = 0

            # for each column chunk in the final matrix
            for j in range(num_whole_chunks):
                # select a section of columns so we get a (n x c) matrix.
                col_mat = geno[:,col_pos_cnt:(col_pos_cnt + cchunk)]

                # select a section of the allele probability column vector
                col_prob = allele_prob[col_pos_cnt:(col_pos_cnt + cchunk),:]

                # calculate the LD in the chunk
                D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                      numpy.matmul(row_prob, col_prob.transpose()))

                # increment col_pos_cnt
                col_pos_cnt += cchunk

                # yield whole thing, otherwise upper triangle matrix.
                yield D

            # calculate a rectangular matrix if one needs to be calculated
            if columns > col_pos_cnt:
                # get everything from col_pos_cnt to end of columns
                col_mat = geno[:,col_pos_cnt:columns]

                # get the last remaining allele probabilities
                col_prob = allele_prob[col_pos_cnt:columns,:]

                # calculate the LD matrix chunk
                D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                      numpy.matmul(row_prob, col_prob.transpose()))

                # yield D. D should never be a square matrix along the diagonal
                # it should always be a rectangular matrix.
                yield D
            # increment the counter by cchunk
            row_pos_cnt += cchunk

        # calculate the last remaining diagonal fi one needs to be calculated
        if columns > row_pos_cnt:
            # get everything from row_pos_cnt to the end
            row_mat = geno[:,row_pos_cnt:columns].transpose()

            # get everything from row_pos_cnt to the end
            row_prob = allele_prob[row_pos_cnt:columns,:]

            # offset the col_pos_cnt variable
            col_pos_cnt = 0

            # for each column chunk in the final matrix
            for j in range(num_whole_chunks):
                # select a section of columns so we get a (n x c) matrix.
                col_mat = geno[:,col_pos_cnt:(col_pos_cnt + cchunk)]

                # select a section of the allele probability column vector
                col_prob = allele_prob[col_pos_cnt:(col_pos_cnt + cchunk),:]

                # calculate the LD in the chunk
                D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                      numpy.matmul(row_prob, col_prob.transpose()))

                # increment col_pos_cnt
                col_pos_cnt += cchunk

                # yield whole thing, otherwise upper triangle matrix.
                yield D

            # calculate a rectangular matrix if one needs to be calculated
            if columns > col_pos_cnt:
                # get everything from col_pos_cnt to end of columns
                col_mat = geno[:,col_pos_cnt:columns]

                # get the last remaining allele probabilities
                col_prob = allele_prob[col_pos_cnt:columns,:]

                # calculate the LD matrix chunk
                D = ((numpy.matmul(row_mat, col_mat) / nseq_f64) -
                      numpy.matmul(row_prob, col_prob.transpose()))

                # yield D. D should never be a square matrix along the diagonal
                # it should always be a rectangular matrix.
                yield D
