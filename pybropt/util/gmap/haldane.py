import numpy

def haldane2r(d):
    """
    Convert genetic map distances to recombination probabilities using the
    Haldane mapping function (Haldane, 1919).

    Parameters
    ----------
    d : numpy.ndarray
        An array of genetic map distances to convert to recombination
        probabilities.

    Returns
    -------
    r : numpy.ndarray
        An array of recombination probabilities between loci as defined by 'd'.
    """
    return 0.5 * (1 - numpy.exp(-2 * d))



def r2haldane(r):
    """
    Convert recombination probabilities between two loci to genetic map
    distances using the Haldane mapping function (Haldane, 1919).

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
