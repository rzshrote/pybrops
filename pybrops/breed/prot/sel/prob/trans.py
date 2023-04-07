"""
Module containing functions for transforming latent space function outputs.
"""

__all__ = [
    "trans_identity",
    "trans_empty"
]

import numpy


def trans_identity(
        vec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Treat all elements in a latent vector as unconstrained objectives.

    Parameters
    ----------
    vec : numpy.ndarray
        An array of shape ``(l,)`` to be transformed.
    kwargs : dict
        Additional keyword arguments. Not used by this function.
    
    Returns
    -------
    out : numpy.ndarray
        An array identical to the input.
    """
    return vec

def trans_empty(
        vec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Return an empty vector.

    Parameters
    ----------
    vec : numpy.ndarray
        An array of shape ``(l,)`` to be transformed.
    kwargs : dict
        Additional keyword arguments. Not used by this function.
    
    Returns
    -------
    out : numpy.ndarray
        An array of shape (0,) with the same dtype as the input.
    """
    return numpy.empty((0,), dtype = vec.dtype)

