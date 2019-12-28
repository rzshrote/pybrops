import numpy

def haldane2r(d):
    """
    Convert genetic distances in Morgans to recombination probabilities using
    the Haldane mapping function (Haldane, 1919).

    This is a bare-bones function and does no parameter checking for errors.

    Parameters
    ==========
    d : numpy.ndarray
        An array of genetic distances in Morgans. This can be an array of any
        shape.

    Returns
    =======
    r : numpy.ndarray
        An array of recombination probabilities. The shape of the array is
        the same shape as that of 'd'.
    """
    return 0.5 * (1 - numpy.exp(-2 * d))

def r2haldane(r):
    """
    Convert recombination probabilities between two loci to genetic map
    distances using the Haldane mapping function (Haldane, 1919).

    This is a bare-bones function and does no parameter checking for errors.

    Parameters
    ----------
    r : numpy.ndarray
        An array of recombination probabilities between two loci.

    Returns
    -------
    d : numpy.ndarray
        An array of genetic map distances as defined by the Haldane mapping
        function.
    """
    return -0.5 * numpy.log(1 - (2 * r))
