import numpy

from .meiosis import meiosis

def wrand(male, female, geno, gmap, lgroup, n,
          pmale = None, pfemale = None,
          interference = None, seed = None):
    """
    Perform a weighted random mating. Males are randomly paired with females
    according to a designated probability of contribution. Note: this means
    that the sample's exact parental contributions are not guaranteed to be
    equivalent to the provided contributions, however the expected contribution
    is guaranteed.

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
    n : int
        Number of progeny to generate.
    pmale : None, numpy.ndarray
        A probability vector of male contributions to the progeny. If None, it
        is assumed that all males are weighted evenly.
    pfemale : None, numpy.ndarray
        A probability vector of the female contributions to the progeny. If
        None, it is assumed that all females are weighted evenly.
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

    # calculate male and female gamete sources
    malesrc = numpy.random.choice(
        male,           # male indices
        size = n,       # select n indices
        replace = True, # with replacement
        p = pmale       # with probability
    )
    femalesrc = numpy.random.choice(
        female,         # female indices
        size = n,       # select n indices
        replace = True, # with replacement
        p = pfemale     # with probability
    )

    # generate male and female gametes.
    malegamete = meiosis(
        geno, gmap, lgroup, malesrc,
        interference = interference
    )
    femalegamete = meiosis(
        geno, gmap, lgroup, femalesrc,
        interference = interference
    )

    # stack our gametes to progeny
    progeny = numpy.stack((malegamete, femalegamete))

    # return our progeny
    return progeny


def erand(male, female, geno, gmap, lgroup,
          emale = None, efemale = None, n = 1,
          interference = None, seed = None):
    """
    Perform an exact contribution random mating. Males and females are assigned
    an exact number of gametes that they contribute to the progeny.

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
    emale : None, numpy.ndarray
        A vector of exact male contributions to the progeny. If None, it
        is assumed that all males contribute one gamete to the offspring.
    efemale : None, numpy.ndarray
        A probability vector of the female contributions to the progeny. If
        None, it is assumed that all females contribute one gamete to the
        offspring.
    n : None, int
        Multiplier for number of progeny to generate.
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

    # if exact male contributions are not provided, assume equal contribution.
    if emale is None:
        emale = numpy.repeat(1, len(male))

    # if exact female contributions are not provided, assume equal contribution.
    if efemale is None:
        efemale = numpy.repeat(1, len(female))

    # generate a list of male and female sources
    malesrc = numpy.repeat(male, n * emale)
    femalesrc = numpy.repeat(female, n * efemale)

    # shuffle the source arrays to make random mating pairs
    numpy.random.shuffle(malesrc)
    numpy.random.shuffle(femalesrc)

    # generate male and female gametes.
    malegamete = meiosis(
        geno, gmap, lgroup, malesrc,
        interference = interference
    )
    femalegamete = meiosis(
        geno, gmap, lgroup, femalesrc,
        interference = interference
    )

    # stack our gametes to progeny
    progeny = numpy.stack((malegamete, femalegamete))

    # return our progeny
    return progeny



def controlled(male, female, geno, gmap, lgroup, n,
               interference = None, seed = None):
    """
    Perform a weighted random mating. Males are randomly paired with females
    according to a designated probability of contribution. Note: this means
    that the sample's exact parental contributions are not guaranteed to be
    equivalent to the provided contributions, however the expected contribution
    is guaranteed.

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

    # calculate male and female gamete sources
    malesrc = numpy.repeat(male, n)
    femalesrc = numpy.repeat(female, n)

    # generate male and female gametes.
    malegamete = meiosis(
        geno, gmap, lgroup, malesrc,
        interference = interference
    )
    femalegamete = meiosis(
        geno, gmap, lgroup, femalesrc,
        interference = interference
    )

    # stack our gametes to progeny
    progeny = numpy.stack((malegamete, femalegamete))

    # return our progeny
    return progeny

def dh(female, geno, gmap, lgroup, n,
       interference = None, seed = None):
    """
    Perform a weighted random mating. Males are randomly paired with females
    according to a designated probability of contribution. Note: this means
    that the sample's exact parental contributions are not guaranteed to be
    equivalent to the provided contributions, however the expected contribution
    is guaranteed.

    Parameters
    ==========
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

    # calculate female gamete source
    femalesrc = numpy.repeat(female, n)

    # generate male and female gametes.
    femalegamete = meiosis(
        geno, gmap, lgroup, femalesrc,
        interference = interference
    )

    # stack two female gametes to get DH inbred progeny
    progeny = numpy.stack((femalegamete, femalegamete))

    # return our progeny
    return progeny
