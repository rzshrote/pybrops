"""
Module containing miscellaneous genetic map subroutines.
"""

__all__ = [
    "cM2d"
]

import numpy

def cM2d(
        cM: numpy.ndarray
    ) -> numpy.ndarray:
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
    d = 0.01 * cM   # convert cM to d
    return d
