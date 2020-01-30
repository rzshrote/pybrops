import numpy
# hack to append into our path the parent directory for this file
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.dirname(os.path.realpath(__file__)))))

from popgen.gmap.const import GMAP_MAP2R_DICT
from popgen.gmap.gmapdist import gmapdist
from util.srange import srange

# population variance for doubled haploid lines
def varA_2way(geno, coeff, gmap, lgroup, t = None,
              mapfn = 'haldane', mem = None):
    """
    Calculate a symmetrical matrix of progeny variance for each pairwise 2-way
    cross according to Osthushenrich et al. (2017). The number of pairwise
    crosses calculated is n(n-1)/2 where n is the number of inbred individuals.

    Parameters
    ==========
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
    coeff : numpy.ndarray
        An array of coefficients for allele effects. The dtype of 'coeff'
        should be a floating point. This array should be single
        dimensional (column,) = (N,) where 'N' represents number of markers.
    gmap : numpy.ndarray
        A 1D array of genetic map positions to convert to recombination
        probabilities. The positions units should be in Morgans.
    lgroup : numpy.ndarray
        A 1D array of linkage group sizes. The sum of the elements should equal
        the length of 'gmap'.
    t : None, int
        Number of generations of random intermating of F1 progeny before doubled
        haploid process. If None, assume 0 generations of random intermating.
        Use None over 0 for increased speed.
        Example:
            F1 -> t = None, (or 0)
            F2 -> t = 1
            random intermating of F2 -> t = 2
    mapfn : None, callable, {'haldane', 'kosambi'}
        Mapping function to use. If None is provided, default to 'haldane'.
        Options:
            callable    : A custom callable Python function.
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
            'haldane'   : Haldane (1919) mapping function
            'kosambi'   : Kosambi (1944) mapping function

    """
    # if 'mem' is None, we are memory greedy and use all we can.
    if mem is None:
        mem = len(gmap)

    # if mapfn is not callable, lookup in available map functions
    if not callable(mapfn):
        mapfn = GMAP_MAP2R_DICT[mapfn] # assume that 'mapfn' is a valid key

    # get the number of individuals
    n_indiv = geno.shape[1]

    # allocate a square matrix for each pairwise variance
    varA = numpy.zeros(
        (n_indiv, n_indiv),
        dtype='float64'
    )

    # linkage group start and stop coordinate variables.
    lst = 0
    lsp = 0

    # for each linkage group
    for g in lgroup:
        # increment 'lsp' by 'g'. This is the stop index (exclusive)
        lsp += g

        # for each computational chunk
        for rst,rsp in zip(range(lst,lsp,mem),srange(lst+mem,lsp,mem)):
            for cst,csp in zip(range(lst,lsp,mem),srange(lst+mem,lsp,mem)):
                # get genetic map distances, calculate recombination probability
                r = mapfn(              # convert distance to probability
                    gmapdist(           # calculate map position distances
                        gmap[rst:rsp],  # gmap for rows
                        gmap[cst:csp]   # gmap for columns
                    )
                )

                # calculate our D matrix; this is specific to our mating scheme
                D = (1.0 - 2.0 * r) / 4.0

                # account for any additional random mating
                if t is not None:
                    D *= (1.0 - r)**t

                # calculate 2x row and col coefficients
                # we need 2x because coefficients are assumed to be per allele.
                rcoeff = 2.0 * coeff[rst:rsp]
                ccoeff = 2.0 * coeff[cst:csp]

                # for each mate pair (excluding selfs)
                for female in range(1,n_indiv): # varA row index
                    for male in range(0,female): # varA col index
                        # calculate our row and col effects vectors
                        # 1) cast to int8 (signed); geno is unsigned.
                        #    The casting should not result in overflow errors
                        #    because we assume geno is a binary matrix.
                        # 2) Take the difference; matrix has values {-1,0,1}
                        # 3) Multiply by 2x row coefficients
                        reffect = (
                            numpy.int8(geno[0,female,rst:rsp]) -
                            numpy.int8(geno[0,male,rst:rsp])
                        ) * rcoeff
                        ceffect = (
                            numpy.int8(geno[0,female,cst:csp]) -
                            numpy.int8(geno[0,male,cst:csp])
                        ) * ccoeff

                        # compute the dot products to get a partial variance
                        varA_part = reffect.dot(D).dot(ceffect)

                        # add this partial variance to the lower triangle
                        varA[female,male] += varA_part

        # set the linkage start for the previous loop to the current stop
        lst = lsp

    # finally now that all variances are calculated for the lower triangle,
    # copy the lower triangle to the upper triangle.
    for female in range(1, n_indiv):
        for male in range(0, female):
            varA[male,female] = varA[female,male]

    # return the symmetrical variance matrix
    return varA
