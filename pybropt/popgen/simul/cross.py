import numpy

import .mate


def f1_dh(male, female, geno, gmap, lgroup, n,
          interference = None, seed = None):
    """
    Simulate double haploids induced from F1 plants.

    Parameters
    ==========
    male : numpy.ndarray
        An array of indices to the 'geno' array indicating the male parents.
    female : numpy.ndarray
        An array of indices to the 'geno' array indicating the female parents.
    geno : numpy.ndarray
        An array of allele states. The dtype of 'geno' should be 'uint8'. Array
        shape should be (depth, row, column) = (M, N, L) where 'M' represents
        number of chromosome phases, 'N' represents number of individuals, 'L'
        represents number of markers. Array format should be the 'C' format.
    gmap : numpy.ndarray
        Genetic map in Morgan units.
    lgroup : numpy.ndarray
        Array of linkage group sizes. The sum of the elements should equal the
        length of 'gmap'.
    n : int, numpy.ndarray
        Number of progeny to generate per cross. If an integer is provided, the
        number of progeny for each cross is equalvalent. If a numpy.ndarray is
        provided, the shape of the array must be the same as male and female.
        Important notes:
            In DH generation, hybrids are first generated from the male and
            female parents according to 'n', then ***1 DH is generated from a
            single hybrid.***
            If you want to induce several DH from one hybrid, write a custom
            function!
    interference : None, int
        Interference term for crossovers.
    seed : None, int
        A seed for initializing the random number generator. If 'seed' is None,
        the internal state of the random number generator is unaltered
        (i.e. the system time is NOT used to seed the RNG engine)

    Returns
    =======
    progeny : numpy.ndarray
        An array of binary allele states in the progeny. Array dimensions are
        (M, n, L) where 'M' represents number of chromosome phases, 'n'
        represents number of progeny specified by the user, and 'L' represents
        number of markers.
    """
    # if a seed is provided, seed the RNG
    if seed is not None:
        numpy.random.seed(seed)

    # make hybrids
    hybrids = mate.controlled(male, female, geno, gmap, lgroup, n)

    # generate DH lines from hybrids
    progeny = mate.dh(
        female = numpy.arange(  # select all hybrids
            hybrids.shape[1]    # number of rows in 'hybrids' matrix
        ),
        geno = hybrids,         # hybrid genotypes
        gmap = gmap,            # genetic map
        lgroup = lgroup,        # linkage group sizes
        n = 1                   # 1 DH per hybrid
    )

    # return the progeny
    return progeny
