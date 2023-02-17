"""
Module containing functions for generating target allele frequencies.
"""

import numpy

__all__ = [
    "target_positive",
    "target_negative",
    "target_stabilizing"
]

def target_positive(u_a: numpy.ndarray) -> numpy.ndarray:
    """
    Generate target allele frequencies, selecting for only beneficial alleles.
    
    Parameters
    ----------
    u_a : numpy.ndarray
        A matrix of additive genomic marker effects of shape ``(p_a,t)``
    
    Returns
    -------
    out : numpy.ndarray
        A matrix of target allele frequencies of shape ``(p_a,t)``
    """
    # desire fixation of beneficial alleles
    return numpy.float64(u_a >= 0.0)

def target_negative(u_a: numpy.ndarray) -> numpy.ndarray:
    """
    Generate target allele frequencies, selecting only deleterious alleles.
    
    Parameters
    ----------
    u_a : numpy.ndarray
        A matrix of additive genomic marker effects of shape ``(p_a,t)``
    
    Returns
    -------
    out : numpy.ndarray
        A matrix of target allele frequencies of shape ``(p_a,t)``
    """
    # desire fixation of deleterious alleles
    return numpy.float64(u_a <= 0.0)

def target_stabilizing(u_a: numpy.ndarray) -> numpy.ndarray:
    """
    Generate target allele frequencies, desiring both alleles in equal proportions.
    
    Parameters
    ----------
    u_a : numpy.ndarray
        A matrix of additive genomic marker effects of shape ``(p_a,t)``
    
    Returns
    -------
    out : numpy.ndarray
        A matrix of target allele frequencies of shape ``(p_a,t)``
    """
    # desire both alleles in equal proportion
    return numpy.full(u_a.shape, 0.5, dtype = "float64")
