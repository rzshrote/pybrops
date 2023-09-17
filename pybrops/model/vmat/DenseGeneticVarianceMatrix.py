"""
Module implementing classes and associated error checking routines for matrices
storing dense genetic variance estimates.
"""

import copy
from typing import Optional
import numpy
import h5py
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_io_h5py import check_group_in_hdf5
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.mat.DenseSquareTaxaTraitMatrix import DenseSquareTaxaTraitMatrix
from pybrops.core.util.h5py import save_dict_to_hdf5
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix

class DenseGeneticVarianceMatrix(DenseSquareTaxaTraitMatrix,GeneticVarianceMatrix):
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
    @epgc.setter
    def epgc(self, value: tuple) -> None:
        """Set data for property epgc."""
        error_readonly("epgc")

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseGeneticVarianceMatrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : DenseGeneticVarianceMatrix
            A shallow copy of the original DenseGeneticVarianceMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseGeneticVarianceMatrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseGeneticVarianceMatrix
            A deep copy of the original DenseGeneticVarianceMatrix.
        """
        return copy.deepcopy(self, memo)

    ################### Matrix File I/O ####################

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> 'DenseGeneticVarianceMatrix':
        """
        Read DenseGeneticVarianceMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which DenseGeneticVarianceMatrix data is stored.
            If None, DenseGeneticVarianceMatrix is read from base HDF5 group.

        Returns
        -------
        out : DenseGeneticVarianceMatrix
            A dense matrix read from file.
        """
        return super(DenseGeneticVarianceMatrix, cls).from_hdf5(
            filename, 
            groupname
        )

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
