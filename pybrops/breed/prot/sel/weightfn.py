"""
Module for building marker weight vectors
"""

import numpy

__all__ = [
    "weight_absolute",
    "weight_one"
]

def weight_absolute(u_a: numpy.ndarray) -> numpy.ndarray:
    """
    Construct marker weights using absolute value of marker effect coefficients
    
    Parameters
    ----------
    u_a : numpy.ndarray
        A matrix of additive genomic marker effects of shape ``(p_a,t)``
    
    Returns
    -------
    out : numpy.ndarray
        A matrix of absolute values of marker effect coefficients of shape ``(p_a,t)``
    """
    return numpy.absolute(u_a)

def weight_one(u_a: numpy.ndarray) -> numpy.ndarray:
    """
    Construct marker weights for even weighing.
    
    Parameters
    ----------
    u_a : numpy.ndarray
        A matrix of additive genomic marker effects of shape ``(p_a,t)``
    
    Returns
    -------
    out : numpy.ndarray
        A matrix of ones of shape ``(p_a,t)``
    """
    return numpy.full(u_a.shape, 1.0, dtype = "float64")
