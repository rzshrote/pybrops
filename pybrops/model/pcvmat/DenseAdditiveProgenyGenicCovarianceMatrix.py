"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense additive genic covariance estimates.
"""

from typing import Optional
import numpy
from pybrops.model.pcvmat.DenseProgenyGenicCovarianceMatrix import DenseProgenyGenicCovarianceMatrix
from pybrops.model.pcvmat.AdditiveProgenyGenicCovarianceMatrix import AdditiveProgenyGenicCovarianceMatrix

class DenseAdditiveProgenyGenicCovarianceMatrix(
        DenseProgenyGenicCovarianceMatrix,
        AdditiveProgenyGenicCovarianceMatrix,
    ):
    """
    A semi-concrete class for dense additive genetic covariance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense additive genetic covariance matrices.

    Methods responsible for estimating genetic variances from genomic models
    remain abstract and must be implemented by inheriting classes.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ):
        """
        Constructor for the concrete class DenseAdditiveProgenyGenicCovarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        trait : numpy.ndarray
            Trait names.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseAdditiveProgenyGenicCovarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    # from_gmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete

    # from_algmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete



################################## Utilities ###################################
def check_is_DenseAdditiveProgenyGenicCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseAdditiveProgenyGenicCovarianceMatrix``.
    Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseAdditiveProgenyGenicCovarianceMatrix):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,DenseAdditiveProgenyGenicCovarianceMatrix.__name__,type(v).__name__))
