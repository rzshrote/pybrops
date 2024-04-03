"""
Module containing transformation subroutines.
"""

import numpy


def trans_ndpt_pseudo_dist(
        ndptmat: numpy.ndarray, 
        objfn_minmax: numpy.ndarray,
        objfn_pseudoweight: numpy.ndarray,
        **kwargs: dict
    ) -> numpy.ndarray:
    """
    Transform a set of non-dominated points by calculating their distances to a
    vector.

    Parameters
    ----------
    ndptmat : numpy.ndarray
        A point coordinate array of shape (npt, nobj) where 'npt' is the number
        of points and 'nobj' is the number of objectives (dimensions).
        This array contains input points for calculating the distance between a
        point to a pseudoweight vector.
        
    objfn_minmax : numpy.ndarray
        An objective function min/max sign indicator array of shape ``(nobj,)``.
        This array contains sign indicator variables for whather an objective is
        minimizing (-1.0) or maximizing (1.0).

    objfn_pseudopreference : numpy.ndarray
        An objective function pseudoweight array of shape ``(nobj,)``

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
    # input assertions
    assert numpy.all(objfn_pseudoweight >= 0.0)
    assert numpy.any(objfn_pseudoweight > 0.0)
    assert objfn_pseudoweight.dot(objfn_pseudoweight) > 0.0
    
    # transform ndptmat to all maximizing functions
    # (npt,nobj) * (nobj,) -> (npt,nobj) * (1,nobj) -> (npt,nobj)
    ndptmat = ndptmat * objfn_minmax

    # subtract column minimums
    # (npt,nobj).min(0) -> (nobj,)
    # (npt,nobj) - (nobj,) -> (npt,nobj) - (1,nobj) -> (npt,nobj)
    ndptmat = ndptmat - ndptmat.min(0)

    # divide by column maximums avoiding division by zero; ndptmat is in range [0,1]
    # (npt,nobj).min(0) -> (nobj,)
    # scalar / (nobj,) -> (nobj,)
    # (nobj,) * (npt,nobj) -> (npt,nobj)
    maximum = ndptmat.max(0)
    mask = (maximum == 0.0)
    maximum[mask] = 1.0
    scale = 1.0 / maximum
    scale[mask] = 0.0
    ndptmat = scale * ndptmat

    # calculate the projection of a point onto the line
    #             P.L
    # proj_L(P) = --- L
    #             L.L
    # where:
    #   P = ndptmat
    #   L = objfn_pseudoweight

    # (nobj,) . (nobj,) -> scalar
    LdotLinv = 1.0 / objfn_pseudoweight.dot(objfn_pseudoweight)
    # (npt,nobj) . (nobj,) -> (npt,)
    PdotL = ndptmat.dot(objfn_pseudoweight)
    # scalar * (npt,) -> (npt,)
    scale = LdotLinv * PdotL
    # (npt,1) * (npt,nobj) -> (npt,nobj)
    projL_P = scale[:,None] * objfn_pseudoweight

    # calculate the orthoprojection of a point onto the line
    # oproj_L(P) = P - proj_L(P)

    # (npt,nobj) - (npt,nobj) -> (npt,nobj)
    oprojL_P = ndptmat - projL_P

    # calculate the euclidean distance of orthoprojection
    # (npt,nobj) -> (npt,)
    dist = numpy.linalg.norm(oprojL_P, axis = 1)

    return dist
