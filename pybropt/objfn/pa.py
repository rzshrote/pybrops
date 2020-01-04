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

################################################################################
# Create a couple of dictionaries to assist the main 'pa' function
################################################################################

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

# population architect selection
def pa(rsel, geno, wcoeff, tfreq, tld, gmap, lgroup, cycles,
       efreq = None, eld = None, mapfn = None, ldfn = None,
       mem = None, mtype = None, dtype = numpy.dtype("float64") ):
    """
    Population Architect (PA) Selection.

    Bare bones objective function. Does minimal checking of arguments.

    Parameters
    ==========
    rsel : numpy.ndarray, list, slice
        A 1D array of indices to use for slicing and summing the matrix by row.
        Each index in 'rsel' represents a single individual's row.
        If 'rsel' == slice(None, None, None) use all individuals.
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
    wcoeff : numpy.ndarray
        An array of coefficients for allele weights. Values in 'wcoeff' should
        be non-negative. The dtype of 'wcoeff' should be either 'float32' or
        'float64'. This array should be of shape (column,) = (N,) where 'N'
        represents number of markers.
    tfreq : numpy.ndarray
        An array of target allele frequencies.
        Example:
            tfreq = numpy.array([0.2, 0.6, 0.7])
    tld : numpy.ndarray
        An array of target LD values. Current measure is in r squared.
        Example:
            tld = numpy.array([1.0, 1.0, 1.0])
    gmap : numpy.ndarray
        Genetic map in Morgan units.
    lgroup : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'gmap'.
    cycles : int, numpy.integer
        Number of breeding cycles left to the deadline.
    efreq : None, numpy.ndarray
        An array of allele frequencies maximally distant from the target
        frequency. If None, function will calculate this (at a speed cost).
        It is better to provide a pre-calculated array if scoring multiple
        selection subsets from the same population and selection criteria.
        Example:
            tfreq   = numpy.array([0.2, 0.6, 0.7])
            efreq   = numpy.array([1.0, 0.0, 0.0])
    eld : None, numpy.ndarray
        An array of LD values maximally distant from the target LD. If None,
        function will calculate this (at a speed cost) assuming a valid LD
        range of [0.0,1.0]. It is better to provide a pre-calculated array if
        scoring multiple selection subsets from the same population and
        selection criteria.
        Example:
            tld     = numpy.array([1.0, 1.0, 0.0])
            eld     = numpy.array([0.0, 0.0, 1.0])
    mapfn : None, callable, {'haldane', 'kosambi'}
        Mapping function to use. If None is provided, default to 'haldane'.
        Options:
            callable  : A custom callable Python function.
                        Must take the following arguments:
                            mapfn(d)
                            Parameters:
                                d : numpy.ndarray
                                    An array of genetic distances in Morgans.
                                    This can be an array of any shape.
                            Returns:
                                r : numpy.ndarray
                                    An array of recombination probabilities.
                                    The shape of the array is the same shape
                                    as that of 'd'.
            'haldane' : Haldane (1919) mapping function
            'kosambi' : Kosambi (1944) mapping function
    ldfn : None, callable, {'r_sq', 'D_prime'}
        Linkage disequilibrium function. Customize this to get custom LD
        measures or crossing schemes.
        Options:
            callable  : A custom callable Python function.
                        Must take the following arguments:
                            ldfn(rgeno, cgeno, rprob, cprob, r, cycles)
                            Parameters: see r_sq function documentation for details
                                rgeno : numpy.ndarray
                                    Row allele states
                                cgeno : numpy.ndarray
                                    Column allele states
                                rprob : numpy.ndarray
                                    Row allele probabilities
                                cprob : numpy.ndarray
                                    Column allele probabilities
                                r : numpy.ndarray
                                    Recombination probabilities
                                cycles : int, numpy.integer
                                    Random mating cycles
            'r_sq'    : R squared linkage disequilibrium metric.
            'D_prime' : D' linkage disequilibrium metric.
    mem : None, {int, numpy.integer}
        Size of memory chunk to compute at a time. This is useful in situations
        with large numbers of markers (>5000). Large matrix calculations are
        exceedingly memory hungry (memory needed = n^2). Set this to a lower
        number to reduce memory consumption during matrix operations.
        Options:
            None:
                Compute everything at once; no memory restrictions.
            int, numpy.int:
                Specify chunks of maximum width 'mem'. Memory is restricted,
                but is proportional to dtype.itemsize * (mem**2)
    mtype : None, {'sym', 'triu', 'tril'}
        Matrix type to return. If None is provided, default to 'sym'
        Options:
            'sym'   : Symmetric matrix
            'triu'  : Upper triangular matrix
            'tril'  : Lower triangular matrix
    dtype : numpy.dtype
        The type of the output score. If dtype is not given, default to
        numpy.dtype('float64').
        Regardless of 'dtype's value, all internal calculations are done with
        double-precision math to reduce floating point errors.

    Returns
    =======
    score : numpy.dtype
        A floating point number representing the Population Architect score for
        the selected subset.
    """
    ############################## Parse 'efreq' ##############################
    # if edge frequency has not been provided, generate it.
    if efreq is None:
        # if tfreq is a numpy.ndarray, generate target matrix
        if isinstance(tfreq, numpy.ndarray):
            efreq = numpy.asarray(
                tfreq < 0.5,        # if tfreq < 0.5 set 1.0; else set 0.0
                dtype=tfreq.dtype   # make identical data types
            )
        # else, we'll assume tfreq is a floating point
        else:
            efreq = 1.0 if tfreq < 0.5 else 0.0

    ############################### Parse 'eld' ###############################
    # if edge linkage disequilibrium has not been provided, generate it.
    if eld is None:
        # if tld is a numpy.ndarray, generate target matrix
        if isinstance(tld, numpy.ndarray):
            eld = numpy.asarray(
                tld < 0.5,          # if tfreq < 0.5 set 1.0; else set 0.0
                dtype=tld.dtype     # make identical data types
            )
        # else, we'll assume tld is a floating point
        else:
            eld = 1.0 if tld < 0.5 else 0.0

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

    # WEIGHT SCORE: calc scaling factor for pairwise matrix
    # divide matrix chunks by this this scaling factor
    divisor = wcoeff.sum(               # take sum of wcoeff array
        dtype=numpy.float64             # force double-precision
    )
    divisor *= divisor                  # square the sum
    if mtype in ["tril", "triu"]:       # if we want a triangle matrix
        divisor += (wcoeff*wcoeff).sum( # add sum of squared elements
            dtype=numpy.float64         # force double-precision
        )
        divisor /= 2.0                  # divide by two to get upper/lower sum

    # ALLELE SCORE, LD SCORE: calculate allele frequencies
    # take sum column and divide by num phases
    pfreq = geno[:,rsel,:].sum((0,1)) / numpy.float64(phases)

    # ALLELE SCORE: calc diff between target and population allele freqs
    diff_tp = tfreq - pfreq

    # ALLELE SCORE: calc diff between target and edge allele freqs
    diff_te = tfreq - efreq

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
            # make exhaustive pairwise product; divide by divisor
            cumprod = (wcoeff[rst:rsp,None] @ wcoeff[None,cst:csp]) / divisor

            #################################################################
            ## Calculate pairwise allele frequency distance scores for all ##
            #################################################################
            # make a mesh for differences between target and population
            mesh_tp = numpy.meshgrid(   # meshgrid returns a list
                diff_tp[rst:rsp],       # n_rows: repeat vector for n_cols (mesh 1)
                diff_tp[cst:csp],       # n_cols: repeat vector for n_rows (mesh 2)
                indexing='ij'           # use 'ij' indexing (faster)
            )

            # make a mesh for differences between edge and population
            mesh_te = numpy.meshgrid(   # meshgrid returns a list
                diff_te[rst:rsp],       # n_rows: repeat vector for n_cols (mesh 1)
                diff_te[cst:csp],       # n_cols: repeat vector for n_rows (mesh 2)
                indexing='ij'           # use 'ij' indexing (faster)
            )

            # calculate distance between target and population allele freq.
            # Gory details:
            #     Recall: diff_tp = tfreq - pfreq        for each locus
            #     Therefore:
            #         mesh_tp[0] = (tfreq_i - pfreq_i)
            #         mesh_tp[1] = (tfreq_j - pfreq_j)
            #     We want to calculate euclidian distance: dist_tp
            #         c = sqrt(a^2 + b^2)
            #         sqrt((tfreq_i - pfreq_i)**2 + (tfreq_j - pfreq_j)**2)
            #     Therefore the calculation is:
            dist_tp = numpy.sqrt(           # sqrt to get euclidian distance
                (mesh_tp[0] * mesh_tp[0]) + # use m*m format (faster)
                (mesh_tp[1] * mesh_tp[1])   # use m*m format (faster)
            )

            # distance between target and edge allele freq is similar:
            # Gory details:
            #     Recall: diff_te = tfreq - efreq
            #     Therefore:
            #         mesh_te[0] = (tfreq_i - efreq_i)
            #         mesh_te[1] = (tfreq_j - efreq_j)
            #     We want to calculate euclidian distance: dist_tp
            #         c = sqrt(a^2 + b^2)
            #         sqrt((tfreq_i - efreq_i)**2 + (tfreq_j - efreq_j)**2)
            #     Therefore the calculation is:
            dist_te = numpy.sqrt(           # sqrt to get euclidian distance
                (mesh_te[0] * mesh_te[0]) + # use m*m format (faster)
                (mesh_te[1] * mesh_te[1])   # use m*m format (faster)
            )

            # divide distance from the target by distance from the edge
            # assuming that the edge is the maximal distance between the
            # target and any hypothetical allele frequency (range [0.0,1.0])
            # this guarantees that the distance metric is in range [0.0,1.0]
            # do (1 - dist_prime) to get an allele score:
            #     allele_score == 0.0  --> maximal distance from target
            #     allele_score == 1.0  --> on target
            cumprod *= (1.0 - (dist_tp / dist_te))

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

            # if: linkage disequilibrium is approximately zero,
            #     the LD component is equal to the allele availability
            # else:
            #     calculate distance scaled distance from target to actual;
            # accumulate product into 'cumprod'
            cumprod *= numpy.where(
                pld < 1e-8,
                allele_avail[rst:rsp,None] @ allele_avail[None,cst:csp],
                1.0 - ((tld - pld) / (tld - eld))
            )

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
