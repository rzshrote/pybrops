"""
Module containing utility functions for determining Pareto efficiency
(non-dominated points).
"""

import numpy

__all__ = [
    "is_pareto_efficient",
]

# based on:
# https://stackoverflow.com/questions/32791911/fast-calculation-of-pareto-front-in-python
def is_pareto_efficient(
        fmat: numpy.ndarray, 
        wt: numpy.ndarray, 
        return_mask: bool = True
    ) -> numpy.ndarray:
    """
    Find pareto-efficient points assuming maximizing functions.

    Parameters
    ----------
    fmat : numpy.ndarray
        A matrix of shape ``(npt, nobj)`` containing fitness values.

        Where:

        - ``npt`` is the number of points
        - ``nobj`` is the number of objectives.
    wt : numpy.ndarray
        Weights to applied to fmat.
    return_mask : bool
        If ``True``, return a mask.

    Returns
    -------
    out : numpy.ndarray
        An array of indices of pareto-efficient points.

        If return_mask is ``True``, this will be an ``(npt,)`` boolean array.
        Otherwise it will be a ``(n_efficient_points, )`` integer array of
        indices.
    """
    fmat = fmat * (wt.flatten()[None,:])    # apply weights
    npt = fmat.shape[0]                     # get number of points
    is_efficient = numpy.arange(npt)        # starting list of efficient points (holds indices)
    pt_ix = 0  # Next index in the is_efficient array to search for
    while pt_ix < len(fmat):
        ndpt_mask = numpy.any(fmat > fmat[pt_ix], axis=1)
        ndpt_mask[pt_ix] = True
        is_efficient = is_efficient[ndpt_mask]  # Remove dominated points
        fmat = fmat[ndpt_mask]
        pt_ix = numpy.sum(ndpt_mask[:pt_ix])+1
    if return_mask:
        is_efficient_mask = numpy.zeros(npt, dtype = bool)
        is_efficient_mask[is_efficient] = True
        return is_efficient_mask
    else:
        return is_efficient
