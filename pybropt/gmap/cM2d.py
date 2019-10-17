def cM2d(cM):
    """
    Convert centiMorgan units to genetic map distance units (equivalent to
    Morgan units).

    Parameters
    ----------
    cM : numpy.ndarray
        An array of centiMorgan map distance values.

    Returns
    -------
    d : numpy.ndarray
        An array of genetic map distance units (Morgans).
    """

    return 0.01*cM
