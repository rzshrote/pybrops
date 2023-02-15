"""
Module containing utility functions for variance estimations.
"""

import numpy

__all__ = ["rprob_filial", "cov_D1s", "cov_D2s", "cov_D1st", "cov_D2st"]

def rprob_filial(r, k):
    """
    Calculate the expected observed recombination rate between two loci at
    filial generation 'k' as specified by Lehermeier (2017).

    Remark:
        Expected observed recombination rate is what would be observed across
        infinite progeny.

    Parameters
    ----------
    r : numpy.ndarray
        Recombination probability matrix.
    k : int, inf
        Selfing filial generation number to derive gametes from.

        +-------------+-------------------------+
        | Example     | Description             |
        +=============+=========================+
        | ``k = 1``   | Derive gametes from F1  |
        +-------------+-------------------------+
        | ``k = 2``   | Derive gametes from F2  |
        +-------------+-------------------------+
        | ``k = 3``   | Derive gametes from F3  |
        +-------------+-------------------------+
        | ``...``     | etc.                    |
        +-------------+-------------------------+
        | ``k = inf`` | Derive gametes from SSD |
        +-------------+-------------------------+

    Returns
    -------
    r_k : numpy.ndarray
        Expected observed recombination rate between two loci at filial
        generation 'k'.
    """
    # multiply matrix by two
    two_r = 2.0 * r

    # calculate first component of r_k
    r_k = two_r / (1.0 + two_r)

    # if k < inf, we do not have SSD and second term is needed
    if k < numpy.inf:
        r_k *= (1.0 - ( (0.5**k) * ((1.0 - two_r)**k) ) )

    # return result
    return r_k

def cov_D1s(r, nself):
    """
    Calculate a D1 matrix.

    This matrix represents the covariance (caused by linkage disequilibrium)
    between genomic loci one recombination event prior to selfing for random
    intermating.

    Parameters
    ----------
    r : numpy.ndarray
        Recombination probability matrix.
    nself : int, inf
        Selfing generation number to derive gametes from.

        +-------------+-------------------------+
        | Example     | Description             |
        +=============+=========================+
        | ``nself = 0``   | Derive gametes from F1  |
        +-------------+-------------------------+
        | ``nself = 1``   | Derive gametes from F2  |
        +-------------+-------------------------+
        | ``nself = 2``   | Derive gametes from F3  |
        +-------------+-------------------------+
        | ``...``     | etc.                    |
        +-------------+-------------------------+
        | ``nself = inf`` | Derive gametes from SSD |
        +-------------+-------------------------+

    Returns
    -------
    D1 : numpy.ndarray
        A D1 matrix of covariance between genomic loci.
    """
    if nself == 0:
        return (1 - 2 * r)
    elif nself > 0:
        return (1.0 - 2.0 * rprob_filial(r, nself+1))
    else:
        raise ValueError("s must be >= 0.")

def cov_D2s(r, nself):
    """
    Calculate a D2 matrix.

    This matrix represents the covariance (caused by linkage disequilibrium)
    between genomic loci two recombination events prior to selfing for random
    intermating.

    Parameters
    ----------
    r : numpy.ndarray
        Recombination probability matrix.
    nself : int, inf
        Selfing generation number to derive gametes from.

        +-------------+-------------------------+
        | Example     | Description             |
        +=============+=========================+
        | ``nself = 0``   | Derive gametes from F1  |
        +-------------+-------------------------+
        | ``nself = 1``   | Derive gametes from F2  |
        +-------------+-------------------------+
        | ``nself = 2``   | Derive gametes from F3  |
        +-------------+-------------------------+
        | ``...``     | etc.                    |
        +-------------+-------------------------+
        | ``nself = inf`` | Derive gametes from SSD |
        +-------------+-------------------------+

    Returns
    -------
    D2 : numpy.ndarray
        A D2 matrix.
    """
    if nself == 0:
        return (1.0 - 2.0 * r)**2
    elif nself > 0:
        four_r = 4.0 * r
        return (1.0 - four_r + (four_r * rprob_filial(r, nself+1)))
    else:
        raise ValueError("s must be >= 0.")

def cov_D1st(r, nself, t):
    """
    Calculate a D1 matrix.

    This matrix represents the covariance between genomic loci caused by
    linkage disequilibrium.

    Parameters
    ----------
    r : numpy.ndarray
        Recombination probability matrix.
    nself : int, inf
        Selfing generation number to derive gametes from.

        +-------------+-------------------------+
        | Example     | Description             |
        +=============+=========================+
        | ``nself = 0``   | Derive gametes from F1  |
        +-------------+-------------------------+
        | ``nself = 1``   | Derive gametes from F2  |
        +-------------+-------------------------+
        | ``nself = 2``   | Derive gametes from F3  |
        +-------------+-------------------------+
        | ``...``     | etc.                    |
        +-------------+-------------------------+
        | ``nself = inf`` | Derive gametes from SSD |
        +-------------+-------------------------+

        Remark:

        - 's' takes priority over 't'. If s > 0, 't' calculations are not
          performed.
    t : int, inf
        Random intermating generation number to derive gametes from.

        +-------------+--------------------------------------------------+
        | Example     | Description                                      |
        +=============+==================================================+
        | ``t = 0``   | Derive gametes from F1                           |
        +-------------+--------------------------------------------------+
        | ``t = 1``   | Derive gametes from (F1 x F1)                    |
        +-------------+--------------------------------------------------+
        | ``t = 2``   | Derive gametes from (F1 x F1) x (F1 x F1)        |
        +-------------+--------------------------------------------------+
        | ``...``     | etc.                                             |
        +-------------+--------------------------------------------------+
        | ``t = inf`` | Derive gametes from unlimited random intermating |
        +-------------+--------------------------------------------------+

        Remark:

        - 's' takes priority over 't'. If s > 0, 't' calculations are not
          performed.

    Returns
    -------
    D_1 : numpy.ndarray
        A D_1 matrix.
    """
    if nself == 0 and t == 0:
        return (1 - 2*r)
    elif nself > 0:
        r_k = rprob_filial(r, nself+1)
        return (1.0 - 2.0*r_k)
    elif t > 0:
        return (1.0 - 2.0*r) * ((1.0 - r)**t)
    else:
        raise ValueError("s and t must be >= 0.")

def cov_D2st(r, nself, t):
    """
    Calculate a D_2 matrix. This matrix represents a linkage disequilibrium
    prototype for two recombination events prior to selfing or random
    intermating.

    Parameters
    ----------
    r : numpy.ndarray
        Recombination probability matrix.
    nself : int, inf
        Selfing generation number to derive gametes from.

        +-------------+-------------------------+
        | Example     | Description             |
        +=============+=========================+
        | ``nself = 0``   | Derive gametes from F1  |
        +-------------+-------------------------+
        | ``nself = 1``   | Derive gametes from F2  |
        +-------------+-------------------------+
        | ``nself = 2``   | Derive gametes from F3  |
        +-------------+-------------------------+
        | ``...``     | etc.                    |
        +-------------+-------------------------+
        | ``nself = inf`` | Derive gametes from SSD |
        +-------------+-------------------------+

        Remark:

        - 's' takes priority over 't'. If s > 0, 't' calculations are not
          performed.
    t : int, inf
        Random intermating generation number to derive gametes from.

        +-------------+--------------------------------------------------+
        | Example     | Description                                      |
        +=============+==================================================+
        | ``t = 0``   | Derive gametes from F1                           |
        +-------------+--------------------------------------------------+
        | ``t = 1``   | Derive gametes from (F1 x F1)                    |
        +-------------+--------------------------------------------------+
        | ``t = 2``   | Derive gametes from (F1 x F1) x (F1 x F1)        |
        +-------------+--------------------------------------------------+
        | ``...``     | etc.                                             |
        +-------------+--------------------------------------------------+
        | ``t = inf`` | Derive gametes from unlimited random intermating |
        +-------------+--------------------------------------------------+

        Remark:

        - 's' takes priority over 't'. If s > 0, 't' calculations are not
          performed.

    Returns
    -------
    D_2 : numpy.ndarray
        A D_2 matrix.
    """
    if nself == 0 and t == 0:
        return (1.0 - 2.0*r)**2
    elif nself > 0:
        r_k = rprob_filial(r, nself+1)
        four_r = 4.0*r
        return ( 1.0 - four_r + (four_r*r_k) )
    elif t > 0:
        return ((1.0 - 2.0*r)**2) * ((1.0 - r)**t)
    else:
        raise ValueError("s and t must be >= 0.")
