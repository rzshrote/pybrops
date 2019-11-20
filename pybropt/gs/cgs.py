# import 3rd party libraries
import numpy
import time
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
import objfn

def cgs(geno,
        coeff,
        k,
        algo = None):
    """
    Conventional Genomic Selection

    Parameters
    ==========
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
        # TODO: performance testing to see which array format is better 'F' or
                'C'.
    coeff : numpy.ndarray
        An array of coefficients for allele effects. The dtype of 'coeff'
        should be either 'float32' or 'float64'. This array should be single
        dimensional (column,) = (N,) where 'N' represents number of markers.
    k : int, numpy.integer, float, numpy.float
        Number of individuals to select OR proportion of individuals to select.
        If a proportion is given, round to the nearest individual.
    algo : None, {‘quicksort’, ‘mergesort’, ‘heapsort’, ‘stable’}
        Optional specification for algorithm to use for sorting GEBVs. Default
        is 'quicksort'. Only worry about this for an extremely large number of
        individuals.

    """
    ########################################
    # Step 0) process arguments
    ########################################

    # check for correct genotype shape
    if len(geno.shape) != 3:
        raise TypeError(
            "'geno' must have 3 dimensions: (phase, indiv, locus)\n"\
            "    geno.shape = %s" %
            (geno.shape,)
        )

    # check for correct coeff shape
    if len(coeff.shape) != 1:
        raise TypeError("'coeff' must have 1 dimension: (locus,)")

    # check if loci lengths are the same
    if geno.shape[2] != coeff.shape[0]:
        raise ValueError(
            "Length of 'coeff' does not match 'geno':\n"\
            "    geno.shape = (%s, %s, locus=%s)\n"\
            "    coeff.shape = (locus=%s,)\n"\
            "'locus' dimensions need to match." %
            (geno.shape[0], geno.shape[1], geno.shape[2], coeff.shape[0])
        )

    # process 'k' and make sure it's the right data type
    if isinstance(k, (int, numpy.integer)):
        k = k if k < geno.shape[1] else geno.shape[1]
    elif isinstance(k, (float, numpy.float)):
        k = numpy.int32(numpy.around(k))
        k = k if k < geno.shape[1] else geno.shape[1]
    else:
        raise TypeError(
            "'k' must be int, numpy.integer, float, or numpy.float\n"\
            "    type(k) = %s" %
            (type(k),)
        )

    # calculate GEBVs
    gebv = objfn.cgs(
        numpy.arange(geno.shape[1]), # select everything
        geno,
        coeff
    )

    # calculate indices for sorted array
    gebv_argsort = numpy.argsort(gebv, kind=algo)

    # return the 'k' highest scoring indices
    return gebv_argsort[-k:].copy()

def cgs_sim(geno,
            d,
            lgroup_size,
            gamete_index,
            interference = None,
            seed = None,
            verbose = True):
    return
