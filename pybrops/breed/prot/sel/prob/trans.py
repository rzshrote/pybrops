"""
Module containing functions for transforming latent space function outputs.
"""

__all__ = [
    "trans_identity",
    "trans_empty",
    "trans_ndpt_to_vec_dist"
]

import numpy


def trans_identity(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Treat all elements in a latent vector as unconstrained objectives.

    Parameters
    ----------
    decnvec : numpy.ndarray
        An array of shape ``(ndecn,)`` containing decision variables 
        corresponding the ``latentvec`` values.
    latentvec : numpy.ndarray
        An array of shape ``(l,)`` to be transformed.
    kwargs : dict
        Additional keyword arguments. Not used by this function.
    
    Returns
    -------
    out : numpy.ndarray
        An array identical to the input.
    """
    return latentvec

def trans_empty(
        decnvec: numpy.ndarray,
        latentvec: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Return an empty vector.

    Parameters
    ----------
    decnvec : numpy.ndarray
        An array of shape ``(ndecn,)`` containing decision variables 
        corresponding the ``latentvec`` values.
    latentvec : numpy.ndarray
        An array of shape ``(l,)`` to be transformed.
    kwargs : dict
        Additional keyword arguments. Not used by this function.
    
    Returns
    -------
    out : numpy.ndarray
        An array of shape (0,) with the same dtype as the input.
    """
    return numpy.empty((0,), dtype = latentvec.dtype)

def trans_ndpt_to_vec_dist(
        mat: numpy.ndarray, 
        obj_wt: numpy.ndarray, 
        vec_wt: numpy.ndarray, 
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Transform a set of non-dominated points by calculating their distances to a
    vector.

    Parameters
    ----------
    mat : numpy.ndarray
        A point coordinate array of shape (npt, nobj) where 'npt' is the number
        of points and 'nobj' is the number of objectives (dimensions).
        This array contains input points for calculating the distance between a
        point to the vector 'wt'.
    obj_wt : numpy.ndarray
        An objective function weight array of shape (nobj,) where 'nobj' is the
        number of objectives (dimensions).
        This array is used to weight the objectives as increasing (1.0) or
        decreasing (-1.0). Objective function weights outside these values have
        undefined behavior.
    vec_wt : numpy.ndarray
        A vector array of shape (nobj,) where 'nobj' is the number of objectives
        (dimensions).
        This array is used as the vector to calculate the distance between each
        point in 'mat'
    kwargs : dict
        Additional keyword arguments. Not used by this function.

    Returns
    -------
    out : numpy.ndarray
        An array of shape (npt,) containing the distance between each point
        to the vector.

    Notes
    -----
    How the transformation is done:

    1. Apply weights to each objective.
    2. Scale points for each objective to the range [0,1]
    3. For each scaled point, calculate a plane that is both perpendicular to
       the input vector and contains the point. Determine the point on the
       plane where the vector intersects the plane.
    4. Calculate the distance (norm) between the two points.

    """
    # transform mat to all maximizing functions
    # (npt,nobj) * (nobj,) -> (npt,nobj) * (1,nobj) -> (npt,nobj)
    mat = mat * vec_wt

    # subtract column minimums
    # (npt,nobj).min(0) -> (nobj,)
    # (npt,nobj) - (nobj,) -> (npt,nobj) - (1,nobj) -> (npt,nobj)
    mat = mat - mat.min(0)

    # divide by column maximums; mat is in range [0,1]
    # (npt,nobj).min(0) -> (nobj,)
    # scalar / (nobj,) -> (nobj,)
    # (nobj,) * (npt,nobj) -> (npt,nobj)
    mat = (1.0 / mat.max(0)) * mat

    # calculate distance between point and line
    # calculate ((v dot p) / (v dot v)) * v
    # where v is the line vector originating from 0
    #       p is the point vector

    # get inverse of (v dot v)
    # (nobj,) dot (nobj,) -> scalar
    vdvinv = 1.0 / obj_wt.dot(obj_wt)

    # get scaling factor (v dot p) / (v dot v)
    # (npt,nobj) dot (nobj,) -> (npt,)
    # (npt,) * scalar -> (npt,)
    scale = mat.dot(obj_wt) * vdvinv

    # use outer product to get points on plane that intersect line
    # (npt,) outer (nobj,) -> (npt,nobj)
    P = numpy.outer(scale, obj_wt)

    # calculate difference between each vector
    # (npt,nobj) - (npt,nobj) -> (npt,nobj)
    diff = mat - P

    # take vector norms of difference to get distances
    # (npt,nobj) -> (npt,)
    d = numpy.linalg.norm(diff, axis = 1)

    return d

