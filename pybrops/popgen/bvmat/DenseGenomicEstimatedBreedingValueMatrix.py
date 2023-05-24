"""
Module implementing matrix routines and associated error checking routines
for dense breeding value matrices estimated from genomic and phenotypic data.
"""

from typing import Any, Optional

import numpy
from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix

# TODO: add standard errors for this class; this could be used for two-stage estimation
class DenseGenomicEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """
    The DenseGenomicEstimatedBreedingValueMatrix class uses a dense matrix to
    represent a Multivariate Breeding Value.

    Notes
    -----
    All elements within a BreedingValueMatrix are mean-centered and scaled to
    unit variance for each trait.

    .. math::
        BV = \\frac{X - \\mu}{\\sigma}

    Where:

    - :math:`BV` is the breeding value.
    - :math:`X` is the phenotype value.
    - :math:`\\mu` is the mean (location) for :math:`X`.
    - :math:`\\sigma` is the standard deviation (scale) for :math:`X`.

    Phenotype values can be reconstituted using:

    .. math::
        X = \\sigma BV + \\mu
    """

    def __init__(
            self, 
            mat: numpy.ndarray, 
            location: numpy.ndarray, 
            scale: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseGenomicEstimatedBreedingValueMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            A numpy.float64 matrix of breeding values of shape ``(n,t)``.
        location : numpy.ndarray
            A numpy.ndarray of shape ``(t,)`` containing breeding value locations.
        scale : numpy.ndarray
            A numpy.ndarray of shape ``(t,)`` containing breeding value scales.
        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DenseGenomicEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGenomicEstimatedBreedingValueMatrix(v: object) -> bool:
    """
    Determine whether an object is a DenseGenomicEstimatedBreedingValueMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseGenomicEstimatedBreedingValueMatrix object instance.
    """
    return isinstance(v, DenseGenomicEstimatedBreedingValueMatrix)

def check_is_DenseGenomicEstimatedBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseGenomicEstimatedBreedingValueMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseGenomicEstimatedBreedingValueMatrix):
        raise TypeError("variable '{0}' must be a DenseGenomicEstimatedBreedingValueMatrix".format(vname))
