"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense genic variance estimates.
"""

from typing import Optional
import numpy
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix

# TODO: implement me
class DenseGenicVarianceMatrix(DenseSquareTaxaMatrix,DenseTraitMatrix,GenicVarianceMatrix):
    """
    UNDER CONSTRUCTION!
    A semi-concrete class for dense genic variance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense genic variance matrices.

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
        Constructor for the concrete class DenseGenicVarianceMatrix.

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

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @DenseSquareTaxaMatrix.square_axes.getter
    def square_axes(self) -> tuple:
        """Get axis indices for axes that are square"""
        return (0,1) # (female, male); same as default in DenseSquareTaxaMatrix

    ######## Expected parental genome contributions ########
    @property
    def epgc(self) -> tuple:
        """Description for property epgc."""
        return (0.5, 0.5)
    @epgc.setter
    def epgc(self, value: tuple) -> None:
        """Set data for property epgc."""
        error_readonly("epgc")
