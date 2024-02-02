"""
Module implementing classes and associated error checking routines for matrices
storing dense genetic covariance estimates.
"""

import copy
from typing import Optional
import numpy
from pybrops.core.mat.DenseSquareTaxaSquareTraitMatrix import DenseSquareTaxaSquareTraitMatrix
from pybrops.model.pcvmat.ProgenyGeneticCovarianceMatrix import ProgenyGeneticCovarianceMatrix


class DenseProgenyGeneticCovarianceMatrix(
        DenseSquareTaxaSquareTraitMatrix,
        ProgenyGeneticCovarianceMatrix,
    ):
    """
    A semi-concrete class for dense genetic covariance matrices

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense genetic covariance matrices.

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
        Constructor for the concrete class DenseProgenyGeneticCovarianceMatrix.

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
        # call DenseSquareTaxaSquareTraitMatrix constructor
        super(DenseProgenyGeneticCovarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @DenseSquareTaxaSquareTraitMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1) # (female, male); same as DenseSquareTaxaSquareTraitMatrix
    
    @DenseSquareTaxaSquareTraitMatrix.square_trait_axes.getter
    def square_trait_axes(self) -> tuple:
        """Axis indices for trait axes that are square."""
        return (2,3) # same as DenseSquareTaxaSquareTraitMatrix

    ######## Expected parental genome contributions ########
    @property
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring from each parent."""
        return (0.5, 0.5)

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseProgenyGeneticCovarianceMatrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : DenseProgenyGeneticCovarianceMatrix
            A shallow copy of the original DenseProgenyGeneticCovarianceMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseProgenyGeneticCovarianceMatrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseProgenyGeneticCovarianceMatrix
            A deep copy of the original DenseProgenyGeneticCovarianceMatrix.
        """
        return copy.deepcopy(self, memo)

    ################### Matrix File I/O ####################

    # to_hdf5               (inherited from DenseSquareTaxaSquareTraitMatrix)

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> 'DenseProgenyGeneticCovarianceMatrix':
        """
        Read DenseProgenyGeneticCovarianceMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which DenseProgenyGeneticCovarianceMatrix data is stored.
            If None, DenseProgenyGeneticCovarianceMatrix is read from base HDF5 group.

        Returns
        -------
        out : DenseProgenyGeneticCovarianceMatrix
            A dense matrix read from file.
        """
        return super(DenseProgenyGeneticCovarianceMatrix, cls).from_hdf5(
            filename, 
            groupname
        )

    # from_gmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete



################################## Utilities ###################################
def check_is_DenseProgenyGeneticCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseProgenyGeneticCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseProgenyGeneticCovarianceMatrix):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,DenseProgenyGeneticCovarianceMatrix.__name__,type(v).__name__))
