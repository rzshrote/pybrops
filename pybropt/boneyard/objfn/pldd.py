import math
import numpy
from itertools import chain
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from util.human2bytes import human2bytes
from util.gmap.haldane import haldane2r
from util.gmap.kosambi import kosambi2r
from util.ld.r_sq_measure import r_sq
from util.ld.D_prime_measure import D_prime
from util.index_factory import ixfty1, ixfty2, ixfty3

MAPFN_DICT = {
    None      : haldane2r,
    'haldane' : haldane2r,
    'kosambi' : kosambi2r
}

LDFN_DICT = {
    None      : r_sq,
    'r_sq'    : r_sq,
    'D_prime' : D_prime
}

MTYPE_IXFTY_DICT = {
    # tuple index 0 = advance row indices: 'start' to 'stop' by 'step'
    # tuple index 1 = advance col indices: 'start' to 'stop' by 'step'
    None    : (ixfty1, ixfty1),
    # tuple index 0 = advance row indices: 'start' to 'stop' by 'step'
    # tuple index 1 = advance col indices: 'start' to 'stop' by 'step'
    "sym"   : (ixfty1, ixfty1),
    # tuple index 0 = advance row indices: 'start' to 'stop' by 'step'
    # tuple index 1 = advance col indices: 'start' to 'rsp' by 'step'
    "tril"  : (ixfty1, ixfty2),
    # tuple index 0 = advance row indices: 'start' to 'stop' by 'step'
    # tuple index 1 = advance col indices: 'rst' to 'stop' by 'step'
    "triu"  : (ixfty1, ixfty3)
}


def pldd(rsel, geno, wcoeff, tfreq, tld, gmap, lgroup, cycles,
         mapfn = 'haldane', ldfn = 'r_sq', mem = None,
         mtype = 'tril', dtype = numpy.dtype('float64')):
    """
    Lower means better.
    """
    ############################## Parse 'mapfn' ##############################
    # if mapfn is not callable(), assume it is a valid key in MAPFN_DICT
    if not callable(mapfn):
        mapfn = MAPFN_DICT[mapfn]

    ############################## Parse 'ldfn' ###############################
    # if ldfn is not callable(), assume it is a valid key in LDFN_DICT
    if not callable(ldfn):
        ldfn = LDFN_DICT[ldfn]

    ############################### Parse 'mem' ###############################
    # if 'mem' is None, set to compute whole matrix
    if mem is None:
        mem = len(gmap)
    # otherwise we assume the provided mem is an int type

    ############################## Parse 'mtype' ##############################
    # search for index generator factory functions
    row_ixfty, col_ixfty = MTYPE_IXFTY_DICT[mtype]
    # if the key is incorrect, dict will raise an error.

    ########################## Compute matrix chunks ##########################

    ###############################################################
    ## create several variables for in matrix calculations below ##

    # ALL: calculate length of gmap
    l = len(gmap)

    # ALL: make variable to store score
    score = numpy.float64(0.0)

    # ALL: calculate number of phases
    phases = geno[:,rsel,:].shape[0] * geno[:,rsel,:].shape[1]

    # LD SCORE: calculate allele frequencies
    # take sum column and divide by num phases
    pfreq = geno[:,rsel,:].sum((0,1)) / numpy.float64(phases)

    # LD SCORE: calculate linkage group stop indices
    lsp_arr = numpy.cumsum(lgroup)

    # LD SCORE: calculate allele availability for the selected population
    allele_avail = numpy.where(
        tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
        pfreq > 0.0,            # then set TRUE if population has >0 allele freq
        numpy.where(            # else
            tfreq > 0.0,        # if target freq > 0.0
            numpy.logical_and(  # then set TRUE if pop freq is between (0.0,1.0)
                pfreq > 0.0,
                pfreq < 1.0
            ),
            pfreq < 1.0         # else set TRUE if pop freq is less than 1.0
        )
    )

    # go through all matrix chunks according to 'mtype' row, col ixfty pattern
    for rst,rsp in row_ixfty(0, l, mem, None, None):    # row loop
        for cst,csp in col_ixfty(0, l, mem, rst, rsp):  # column loop
            # NOTE: cumprod is a variable that accumulates fractional score
            # components for each marker; we use cumprod *= matrix instead
            # of making several matrices to have element-wise products
            # calculated to reduce memory usage by having garbage collector
            # eliminate intermediate matrices

            ################################################################
            ## Calculate pairwise marker weights for the weight component ##
            ################################################################
            # make pairwise product; store into cumulative product matrix
            cumprod = wcoeff[rst:rsp,None] @ wcoeff[None,cst:csp]

            ####################################################################
            ## Calculate linkage disequilibrium distance scores for all pairs ##
            ####################################################################
            # recombination probabilities
            r = numpy.empty(            # allocate empty matrix
                (rsp-rst, csp-cst),     # (rows, columns) dimensions
                dtype=numpy.float64     # set data type
            )

            # variables to track columns we haven't filled yet
            lcols = csp  # cols to the left not filled (index)
            rcols = cst  # cols to the right not filled (index)

            # Algorithm:
            #     For each linkage group, calculate the overlap with the
            #     current computational chunk. If there is an overlap,
            #     compute the overlap. Otherwise, fill with 0.5.

            # start coordinate variable for the linkage groups.
            lst = 0

            # cycle through all the linkage group sectors trying to identify
            # regions of overlap between the computational chunk and the
            # linkage group block
            for lsp in lsp_arr:
                # calculate an overlap region defining a box: (x1,y1,x2,y2)
                # where x1 < x2, y1 < y2
                x1 = max(lst, cst)   # x1; matrix start column
                y1 = max(lst, rst)   # y1; matrix start row
                x2 = min(lsp, csp)   # x2; matrix stop column
                y2 = min(lsp, rsp)   # y2; matrix start row

                # if we've found a linkage group that overlaps
                if (x2 >= x1) and (y2 >= y1):
                    # get the cols overlap (x range) as a slice
                    cols = slice(x1-cst, x2-cst)

                    # make distance mesh for the overlap region
                    # first argument specifies number of columns
                    # second argument specifies number of rows
                    dmesh = numpy.meshgrid(gmap[y1:y2], gmap[x1:x2], indexing='ij')

                    # fill the recombination frequencies in the matrix
                    #   rst-y1:rst-y2     # rows (y)
                    #   cols              # columns (x)
                    r[y1-rst:y2-rst,cols] = mapfn(
                        numpy.absolute(dmesh[0] - dmesh[1])
                    )

                    # fill rows above the linkage block with 0.5
                    r[:y1-rst,cols].fill(0.5)

                    # fill rows below the linkage block with 0.5
                    r[y2-rst:,cols].fill(0.5)

                    # alter column tracking variables
                    if x1 < lcols:  # if x1 is before lcols index
                        lcols = x1  # set the index to x1
                    if x2 > rcols:  # if x2 is after rcols index
                        rcols = x2  # set the index to x2

                # advance the linkage group start index to the current stop
                lst = lsp

            # fill the remaining untouched columns with 0.5
            if lcols > rcols:             # if lcols and rcols overlap
                r.fill(0.5)               # fill the whole thing in one go
            else:                         # else fill each side
                r[:,:lcols-cst].fill(0.5) # fill left
                r[:,rcols-cst:].fill(0.5) # fill right

            # calculate population LD
            # NOTE: do not pass view of matrix view: it creates a copy
            pld = ldfn(
                geno[:,rsel,rst:rsp],   # row markers
                geno[:,rsel,cst:csp],   # col markers
                pfreq[rst:rsp],         # row frequencies
                pfreq[cst:csp],         # col frequencies
                r,                      # recombination probabilities
                cycles                  # number of random matings after
            )

            # calculate the difference in LD between the target and population
            diff_ld = numpy.absolute(tld - pld)

            # extract nan indices for everything
            nan_mask = numpy.isnan(diff_ld)

            # calculate pairwise availabilities
            pair_avail = allele_avail[rst:rsp,None] @ allele_avail[None,cst:csp]

            # if nan, set to its pairwise allele availability (1=yes, 0=no)
            diff_ld[nan_mask] = pair_avail[nan_mask]

            # multiply diff_ld matrix to cumulative product
            cumprod *= diff_ld

            # add triangle modifications if necessary
            if (mtype == 'tril') and (rst == cst):
                cumprod *= numpy.tri(       # multiply by mask
                    *cumprod.shape[-2:],    # get first two dimensions
                    k=0,                    # offset to get lower triangle
                    dtype=cumprod.dtype     # use same dtype; marginal speed up
                )
            elif (mtype == 'triu') and (rst == cst):
                cumprod *= (                    # multiply by mask
                    cumprod.dtype.type(1.0) -   # invert triangle mask for triu
                    numpy.tri(                  # make inverse of triu mask
                        *cumprod.shape[-2:],    # get first two dimensions
                        k=-1,                   # offset to get lower triangle
                        dtype=cumprod.dtype     # use same data type
                    )
                )

            # finally, accumulate the matrix and add to score
            score += cumprod.sum()

    # return the score casted as the output data type
    return dtype.type(score)
