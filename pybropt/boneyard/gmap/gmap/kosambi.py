import numpy

def kosambi2r(d):
    """
    Convert genetic map distances to recombination probabilities using the
    Kosambi mapping function (Kosambi, 1944).

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
    return 0.5 * numpy.tanh(2 * d)



def r2kosambi(r):
    """
    Convert recombination probabilities between two loci to genetic map
    distances using the Kosambi mapping function (Kosambi, 1944).

    Parameters
    ----------
    r : numpy.ndarray
        An array of recombination probabilities between two loci.

    Returns
    -------
    d : numpy.ndarray
        An array of genetic map distances as defined by the Kosambi mapping
        function.
    """
    return numpy.log(1 + (2 * r)) / (4 - (8 * r))
