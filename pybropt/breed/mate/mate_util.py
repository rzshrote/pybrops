def mat_meiosis(geno, sel, xoprob, seed = None):
    """
    Perform meiosis on matrix inputs.

    Parameters
    ----------
    geno : numpy.ndarray
    sel : numpy.ndarray
    xoprob : numpy.ndarray
    seed : int

    Returns
    -------
    gamete : numpy.ndarray
    """
    # if there is a RNG seed, seed the RNG
    pybropt.util.cond_seed_rng(seed)

    # calculate shape of the retured gamete matrix
    # number of rows is the number of elements in sel
    gshape = (len(xoprob), len(sel))

    # generate random numbers to determine crossover points:
    # generate in interval [0,1)
    rnd = numpy.random.uniform(0, 1, gshape)

    # allocate gamete array
    gamete = numpy.empty(gshape, dtype = geno.dtype)

    for i,s in enumerate(sel):
        # calculate locations where crossover occurs
        xoix = numpy.flatnonzero(rnd[i] < xoprob)

        # get starting phase
        phase = 0

        # starting index for copy
        stix = 0

        for spix in xoix:
            # copy from start index (inclusive) to stop index (exclusive)
            gamete[i,stix:spix] = geno[phase,s,stix:spix]

            # move the start index in the next iteration to the current stop index
            stix = spix

            # alternate between phases
            phase = 1 - phase

        # finally copy the last remaining chromosome segments
        gamete[i,stix:] = geno[phase,s,stix:]

    return gamete

def mat_mate(fgeno, mgeno, fsel, msel, xoprob, seed = None):
    """
    Perform mating on matrix inputs.

    Parameters
    ----------
    geno : numpy.ndarray
    fsel : numpy.ndarray
    msel : numpy.ndarray
    xoprob : numpy.ndarray
    seed : int

    Returns
    -------
    progeny : numpy.ndarray
    """
    # if there is a RNG seed, seed the RNG
    pybropt.util.cond_seed_rng(seed)

    # generate gametes
    fgamete = mat_meiosis(fgeno, fsel, xoprob)
    mgamete = mat_meiosis(mgeno, msel, xoprob)

    # generate offspring genotypes by stacking matrices to make 3d matrix
    progeny = numpy.stack([fgamete, mgamete])

    return progeny

def mat_dh(geno, sel, xoprob, seed = None):
    """
    Perform doubled haploid production on matrix inputs.

    Parameters
    ----------
    geno : numpy.ndarray
    sel : numpy.ndarray
    xoprob : numpy.ndarray
    seed : int

    Returns
    -------
    progeny : numpy.ndarray
    """
    # if there is a RNG seed, seed the RNG
    pybropt.util.cond_seed_rng(seed)

    # generate gametes
    gamete = mat_meiosis(geno, sel, xoprob)

    # generate offspring genotypes by stacking matrices to make 3d matrix
    progeny = numpy.stack([gamete, gamete])

    return progeny