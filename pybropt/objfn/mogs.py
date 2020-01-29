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

def mogs(rsel, geno, wcoeff, tfreq, dcoeff, dtype = numpy.dtype('float64')):
    """
    Multi-objective genomic selection objective function.
        The goal is to minimize this function. Lower is better.
        This is a bare bones function. Minimal error checking is done.

    Given a 2D weight vector 'dcoeff', calculate the Euclidian distance from the
    origin according to:
        dist = dot( dcoeff, F(x) )
        Where:
            F(x) is a vector of objective functions:
                F(x) = < f_PAU(x), f_PAFD(x) >

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
    """
    ###############################################################
    ## create several variables for in matrix calculations below ##
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
    F_x = numpy.zeros(2, dtype='float64')
    ###############################################################

    ############################# Compute f_PAU(x) #############################
    F_x[0] = numpy.logical_not(allele_avail).dot(wcoeff)

    ############################# Compute f_PAU(x) #############################
    F_x[1] = numpy.absolute(tfreq - pfreq).dot(wcoeff)

    ###################### Compute F(x) and return score #######################
    # multiply by dimension weights
    F_x *= dcoeff

    # calculate distance from center
    score = numpy.sqrt(F_x.dot(F_x))

    # return the score casted as the output data type
    return dtype.type(score)


def mogs_max(rsel, geno, wcoeff, tfreq, dcoeff, dtype = numpy.dtype('float64')):
    """
    Negate Allele Reaper (AR) selection function 'pa_v2' so that the
    objective becomes to maximize the function instead.
    """
    return -1 * mogs(**locals())
