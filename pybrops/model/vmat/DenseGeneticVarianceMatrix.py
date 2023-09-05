"""
Module implementing classes and associated error checking routines for matrices
storing dense genetic variance estimates.
"""

from typing import Optional

import numpy
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix

class DenseGeneticVarianceMatrix(DenseSquareTaxaMatrix,DenseTraitMatrix,GeneticVarianceMatrix):
    """
    A semi-concrete class for dense genetic variance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense genetic variance matrices.

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
        ) -> None:
        """
        Constructor for the concrete class DenseGeneticVarianceMatrix.

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


    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    # from_gmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete



################################## Utilities ###################################
def check_is_DenseGeneticVarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a DenseGeneticVarianceMatrix".format(vname))
