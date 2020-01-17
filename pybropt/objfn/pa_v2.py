import math
import numpy
from itertools import chain
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
from util.human2bytes import human2bytes
from util.gmap.haldane import haldane2r
from util.gmap.kosambi import kosambi2r
from util.ld.r_sq_measure import r_sq_imp
from util.ld.D_prime_measure import D_prime
from util.index_factory import ixfty1, ixfty2, ixfty3

MAPFN_DICT = {
    None      : haldane2r,
    'haldane' : haldane2r,
    'kosambi' : kosambi2r
}

LDFN_DICT = {
    None      : r_sq_imp,
    'r_sq'    : r_sq_imp#,
    #'D_prime' : D_prime
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

def pa_v2(rsel, geno, wcoeff, tfreq, tld, gmap, lgroup, generations, dcoeff,
          mapfn = 'haldane', ldfn = 'r_sq', mtype = 'tril', mem = None, 
          dtype = numpy.dtype('float64')):
    """
    Population Architect v2 (PA v2) objective function.
        The goal is to minimize this function. Lower is better.
        This is a bare bones function. Minimal error checking is done.

    Given a 3D weight vector 'dcoeff', calculate the Euclidian distance from the
    origin according to:
        dist = dot( dcoeff, F(x) )
        Where:
            F(x) is a vector of objective functions:
                F(x) = < f_PAU(x), f_PAFD(x), f_PLDD(x) >

    f_PAU(x):

    Given the provided genotype matrix 'geno' and row selections from it 'rsel',
    calculate the selection allele freq. From the selection allele frequencies
    and the target allele frequencies, determine if the target frequencies
    cannot be attained after unlimited generations and selection rounds.
    Multiply this vector by a weight coefficients vector 'wcoeff'.

    f_PAFD(x):

    Given a genotype matrix, a target allele frequency vector, and a vector of
    weights, calculate the distance between the selection frequency and the
    target frequency.

    f_PLDD(x):

    Calculate a lower triangle linkage disequilibrium matrix for all provided
    loci. Calculate the distance between a target LD and the population LD.
    Multiply the difference by pairwise weight and take the sum across the
    matrix.
                PLDD = w^T @ tril(abs(T - P)) @ w
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
    # otherwise we assume the provided mem is an int type
    if mem is None:
        mem = len(gmap)

    ############################## Parse 'mtype' ##############################
    # search for index generator factory functions
    # if the key is incorrect, dict will raise an error.
    row_ixfty, col_ixfty = MTYPE_IXFTY_DICT[mtype]

    ###############################################################
    ## create several variables for in matrix calculations below ##
    # process genetic map components
    l = len(gmap)                   # get length of gmap
    lsp_arr = numpy.cumsum(lgroup)  # calculate linkage group stop indices

    # make variable to store cumulative score
    score = numpy.float64(0.0)

    # calculate number of phases
    phases = geno[:,rsel,:].shape[0] * geno[:,rsel,:].shape[1]

    # calculate allele frequencies
    pfreq = geno[:,rsel,:].sum((0,1)) / numpy.float64(phases)

    # calculate allele availability for the selected population
    pfreq_gt_0 = (pfreq > 0.0)  # make mask for pfreq > 0
    pfreq_lt_1 = (pfreq < 1.0)  # make mask for pfreq < 1
    allele_avail = numpy.where( # assess allele availability
        tfreq >= 1.0,           # if target freq >= 1.0 (should always be 1.0)
        pfreq_gt_0,             # then set TRUE if population has >0 allele freq
        numpy.where(            # else
            tfreq > 0.0,        # if target freq > 0.0
            numpy.logical_and(  # then set TRUE if pop freq is between (0.0,1.0)
                pfreq_gt_0,
                pfreq_lt_1
            ),
            pfreq_lt_1          # else set TRUE if pop freq is less than 1.0
        )
    )

    # make a vector of doubles to store scores for F(x)
    F_x = numpy.zeros(3, dtype='float64')
    ###############################################################

    ############################# Compute f_PAU(x) #############################
    F_x[0] = numpy.logical_not(allele_avail).dot(wcoeff)

    ############################# Compute f_PAU(x) #############################
    F_x[1] = numpy.absolute(tfreq - pfreq).dot(wcoeff)

    ####################### Compute f_PLDD(x) in chunks ########################
    # go through all matrix chunks according to 'mtype' row, col ixfty pattern
    for rst,rsp in row_ixfty(0, l, mem, None, None):    # row loop
        for cst,csp in col_ixfty(0, l, mem, rst, rsp):  # column loop
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
                # calculate an overlap region defining a box: (c1,r1,c2,r2)
                # where c1 < c2, r1 < r2
                c1 = max(lst, cst)   # c1; matrix start column
                r1 = max(lst, rst)   # r1; matrix start row
                c2 = min(lsp, csp)   # c2; matrix stop column
                r2 = min(lsp, rsp)   # r2; matrix start row

                # if we've found a linkage group that overlaps
                if (c2 >= c1) and (r2 >= r1):
                    # get the cols overlap (x range) as a slice
                    cols = slice(c1-cst, c2-cst)

                    # make distance mesh for the overlap region
                    dmesh = numpy.meshgrid(
                        gmap[r1:r2],    # number of rows
                        gmap[c1:c2],    # number of columns
                        indexing='ij'   # use ij indexing for speed
                    )

                    # fill the recombination frequencies in the matrix
                    #   rst-r1:rst-r2     # rows (y)
                    #   cols              # columns (x)
                    r[r1-rst:r2-rst,cols] = mapfn(
                        numpy.absolute(dmesh[0] - dmesh[1])
                    )

                    # fill rows above the linkage block with 0.5
                    r[:r1-rst,cols].fill(0.5)

                    # fill rows below the linkage block with 0.5
                    r[r2-rst:,cols].fill(0.5)

                    # alter column tracking variables
                    if c1 < lcols:  # if c1 is before lcols index
                        lcols = c1  # set the index to c1
                    if c2 > rcols:  # if c2 is after rcols index
                        rcols = c2  # set the index to c2

                # advance the linkage group start index to the current stop
                lst = lsp

            # fill the remaining untouched columns with 0.5
            if lcols > rcols:             # if lcols and rcols overlap
                r.fill(0.5)               # fill the whole thing in one go
            else:                         # else fill each side
                r[:,:lcols-cst].fill(0.5) # fill left
                r[:,rcols-cst:].fill(0.5) # fill right

            # calculate population LD ***with imputations***
            # NOTE: do not pass view of matrix view: it creates a copy
            pld = ldfn(
                geno[:,rsel,rst:rsp],   # row markers
                geno[:,rsel,cst:csp],   # col markers
                pfreq[rst:rsp],         # row frequencies
                pfreq[cst:csp],         # col frequencies
                allele_avail[rst:rsp]   # row allele availabilities
                allele_avail[cst:csp]   # column allele availabilities
                r,                      # recombination probabilities
                generations             # random matings & meiosis after
            )

            # calculate LD difference between target and population
            diff_ld = numpy.absolute(
                tld[rst:rsp,cst:csp] - pld
            ) if isinstance(
                tld, numpy.ndarray
            ) else numpy.absolute(tld - pld)

            # add triangle modifications if necessary
            if (rst == cst) and (mtype == 'tril'):
                diff_ld *= numpy.tri(       # multiply by mask
                    *diff_ld.shape[-2:],    # get first two dimensions
                    k=0,                    # offset to get lower triangle
                    dtype=diff_ld.dtype     # use same dtype; marginal speed up
                )
            elif (rst == cst) and (mtype == 'triu'):
                diff_ld *= (                    # multiply by mask
                    diff_ld.dtype.type(1.0) -   # invert triangle mask for triu
                    numpy.tri(                  # make inverse of triu mask
                        *diff_ld.shape[-2:],    # get first two dimensions
                        k=-1,                   # offset to get lower triangle
                        dtype=diff_ld.dtype     # use same data type
                    )
                )
            ####################################################################
            # finally, take w^T @ tril(abs(T - P)) @ w and add to score
            ####################################################################
            F_x[2] += wcoeff[rst:rsp].dot(diff_ld).dot(wcoeff[cst:csp])

    ###################### Compute F(x) and return score #######################
    # multiply by dimension weights
    F_x *= dcoeff

    # calculate distance from center
    score = numpy.sqrt(F_x.dot(F_x))

    # return the score casted as the output data type
    return dtype.type(score)
