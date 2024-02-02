"""
Module implementing classes and associated error checking routines for matrices
storing dense genic covariance estimates.
"""

import copy
from typing import Optional
import numpy
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.model.pcvmat.ProgenyGenicCovarianceMatrix import ProgenyGenicCovarianceMatrix

class DenseProgenyGenicCovarianceMatrix(
        DenseSquareTaxaTraitMatrix,
        ProgenyGenicCovarianceMatrix,
    ):
    """
    A semi-concrete class for dense genic covariance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense genic covariance matrices.

    Methods responsible for estimating genic variances from genomic models
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
        ) -> None:
        """
        Constructor for the concrete class DenseProgenyGenicCovarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        # since this is multiple inheritance, do not use parental constructors
        self.mat = mat
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.trait = trait
        # set taxa metadata to None
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @DenseSquareTaxaTraitMatrix.square_axes.getter
    def square_axes(self) -> tuple:
        """Get axis indices for axes that are square"""
        return (0,1) # (female, male); same as default in DenseSquareTaxaTraitMatrix

    ######## Expected parental genome contributions ########
    @property
    def epgc(self) -> tuple:
        """Description for property epgc."""
        return (0.5, 0.5)

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseProgenyGenicCovarianceMatrix':
        """
        Make a shallow copy of the DenseProgenyGenicCovarianceMatrix.

        Returns
        -------
        out : DenseProgenyGenicCovarianceMatrix
            A shallow copy of the original DenseProgenyGenicCovarianceMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseProgenyGenicCovarianceMatrix':
        """
        Make a deep copy of the DenseProgenyGenicCovarianceMatrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseProgenyGenicCovarianceMatrix
            A deep copy of the original DenseProgenyGenicCovarianceMatrix.
        """
        return copy.deepcopy(self, memo)

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################

    # from_gmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete



################################## Utilities ###################################
def check_is_DenseProgenyGenicCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseProgenyGenicCovarianceMatrix``.
    Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseProgenyGenicCovarianceMatrix):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,DenseProgenyGenicCovarianceMatrix.__name__,type(v).__name__))
