def pop_ld_cross(pRec, mpAB, fpAB, mfreq, ffreq):
    """
    Function from the boneyard to calculate a population x population cross.
    Calculate LD resulting from a cross between two populations.

    Parameters
    ----------
    pRec : numpy.ndarray
        Probability of recombination.
    mpAB : numpy.ndarray
        Male P(AB) matrix.
    fpAB : numpy.ndarray
        Female P(AB) matrix.
    mfreq : numpy.ndarray
        Male allele frequency matrix.
    ffreq : numpy.ndarray
        Female allele frequency matrix.

    Returns
    -------
    pop_ld : numpy.ndarray
    """
    coupling = 0.5 * (1.0 - pRec) * (mpAB + fpAB)
    repulsion = (0.5 * pRec * ((mfreq[:,None] * ffreq) + (ffreq[:,None] * mfreq)))

    pop_ld = coupling + repulsion

    return pop_ld
