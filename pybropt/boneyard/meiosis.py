import numpy
import time

def meiosis(geno,
            gmap,
            lgroup,
            gamete_index,
            interference = None,
            seed = None,
            verbose = True):
    """
    Simulate meiosis. Generate gametes from a genotype matrix.
    NOTE: only handles diploids

    Parameters
    ==========
    geno : numpy.ndarray
    gmap : numpy.ndarray
    lgroup : numpy.ndarray
    gamete_index : numpy.ndarray
    interference : None or int
    seed : None or int
    verbose : boolean
    """

    # if there is no seed, use time to seed rng
    if seed == None:
        seed = numpy.uint32(time.time())

    # seed the rng
    numpy.random.seed(seed)

    if verbose:
        print("Meiosis engine state:")
        print("Genotype array at:", id(geno))
        print("Genotype shape =", geno.shape)
        print("Linkage group sizes =", lgroup)
        print("Gamete index =", gamete_index)
        print("Interference =", interference)
        print("numpy.random seed =", seed)

    # make an empty matrix to contain gametes
    gamete = numpy.empty(
        (len(gamete_index),geno.shape[2]),
        dtype=geno.dtype
    )

    # calculate stop indices for linkage group maps
    dsp = numpy.cumsum(lgroup)

    # calculate start indices for linkage group maps
    dst = dsp - lgroup[0]

    ##################################################
    # calculate lambdas for the Poisson distribution #
    ##################################################

    # allocate memory for an array of lambda values
    # these lambda values will be used to determine how many recombination
    # points to generate from the Poisson distribution
    lambdas = numpy.empty(
        len(lgroup),
        dtype=gmap.dtype
    )

    # for start, stop indices subtract gmap[start] from gmap[stop]
    # this calculates the lambda for each linkage group
    for i,(st,sp) in enumerate(zip(dst, dsp)):
        lambdas[i] = gmap[sp-1] - gmap[st]

    # generate Poisson samples to determine the number of crossover events
    # each row corresponds to a gamete_index; each column, a linkage group
    chiasmata = numpy.random.poisson(
        lambdas,
        (len(gamete_index), len(lgroup))
    )

    # generate phase decisions for each gamete and linkage group
    # 0 will corespond to slice 0; 1 will corespond to slice 1.
    # matrix rows correspond to gametes; each column, a linkage group
    phase = numpy.random.binomial(
        1,
        0.5,
        (len(gamete_index),len(lgroup))
    )

    # for each gamete
    for i in range(len(gamete_index)):
        # for each linkage group
        for j,(st,sp) in enumerate(zip(dst,dsp)):
            # TODO: support for interference, what is implemented is no interference

            # generate a number of chiasmata points between map start, stop
            xo_pts = numpy.sort(
                numpy.random.uniform(
                    gmap[st],
                    gmap[sp-1],
                    chiasmata[i,j]
                )
            )

            # for each value in the xo_pts array
            # if no chiasma were generated, this loop skips
            for chiasma in xo_pts:
                # first do binary search for map chiasmata index
                cindex = numpy.searchsorted(
                    gmap[st:sp],
                    chiasma
                )

                # copy a chosen phase fragment from individual
                gamete[i, st:st+cindex] = \
                    geno[phase[i,j], gamete_index[i], st:st+cindex]

                # increment the st by cindex
                st += cindex

                # alternate the phase for the next loop iteration
                phase[i,j] = 1 - phase[i,j]

            # finish up by copying all remaining sites to the gamete.
            # if there are no chiasma, this is the entire chromosome.
            gamete[i, st:sp] = \
                geno[phase[i,j], gamete_index[i], st:sp]

    # return 2d array of gametes
    return gamete