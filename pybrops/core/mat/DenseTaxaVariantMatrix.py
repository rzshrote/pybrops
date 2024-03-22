"""
Module defining implementing dense matrices with taxa and variant metadata and
associated error checking routines.
"""

import copy
from pathlib import Path
import numpy
from typing import Sequence
from typing import Union
from typing import Optional
from numpy.typing import ArrayLike
import h5py

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable
from pybrops.core.error.error_value_numpy import check_ndarray_ndim_gteq
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.mat.DenseTaxaMatrix import DenseTaxaMatrix
from pybrops.core.mat.DenseVariantMatrix import DenseVariantMatrix
from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybrops.core.util.h5py import h5py_File_read_ndarray
from pybrops.core.util.h5py import h5py_File_read_ndarray_utf8
from pybrops.core.util.h5py import h5py_File_write_dict

class DenseTaxaVariantMatrix(
        DenseTaxaMatrix,
        DenseVariantMatrix,
        TaxaVariantMatrix
    ):
    """
    A concrete class for dense matrices with taxa and variant metadata.

    The purpose of this concrete class is to merge the following implementations 
    and interfaces:

        1. DenseTaxaMatrix (implementation)
        2. DenseVariantMatrix (implementation)
        3. TaxaVariantMatrix (interface)
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None,
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None,
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None,
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for creating dense matrices with taxa and variant metadata.

        Parameters
        ----------
        mat : numpy.ndarray
            A numpy.ndarray used to construct the object.
        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        vrnt_chrgrp : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            group labels. If ``None``, do not store any variant chromosome group
            label information.
        vrnt_phypos : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            physical positions. If ``None``, do not store any variant chromosome
            physical position information.
        vrnt_name : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant names.
            If ``None``, do not store any variant names.
        vrnt_genpos : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            genetic positions. If ``None``, do not store any variant chromosome
            genetic position information.
        vrnt_xoprob : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant crossover
            probabilities. If ``None``, do not store any variant crossover
            probabilities.
        vrnt_hapgrp : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant haplotype
            group labels. If ``None``, do not store any variant haplotype group
            label information.
        vrnt_hapalt : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant alternative
            alleles. If ``None``, do not store any variant alternative allele
            information.
        vrnt_hapref : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant reference
            alleles. If ``None``, do not store any variant reference allele
            information.
        vrnt_mask : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing a variant mask.
            If ``None``, do not store any variant mask information.
        kwargs : dict
            Additional keyword arguments.
        """
        # cannot use super constructor since DenseTaxaMatrix and DenseVariantMatrix
        # constructors have overlapping arguments

        # set data
        self.mat = mat
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.vrnt_chrgrp = vrnt_chrgrp
        self.vrnt_phypos = vrnt_phypos
        self.vrnt_name = vrnt_name
        self.vrnt_genpos = vrnt_genpos
        self.vrnt_xoprob = vrnt_xoprob
        self.vrnt_hapgrp = vrnt_hapgrp
        self.vrnt_hapalt = vrnt_hapalt
        self.vrnt_hapref = vrnt_hapref
        self.vrnt_mask = vrnt_mask

        # set taxa metadata
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

        # set variant metadata
        self.vrnt_chrgrp_name = None
        self.vrnt_chrgrp_stix = None
        self.vrnt_chrgrp_spix = None
        self.vrnt_chrgrp_len = None

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseTaxaVariantMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A copy of the DenseTaxaVariantMatrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            vrnt_chrgrp = copy.copy(self.vrnt_chrgrp),
            vrnt_phypos = copy.copy(self.vrnt_phypos),
            vrnt_name = copy.copy(self.vrnt_name),
            vrnt_genpos = copy.copy(self.vrnt_genpos),
            vrnt_xoprob = copy.copy(self.vrnt_xoprob),
            vrnt_hapgrp = copy.copy(self.vrnt_hapgrp),
            vrnt_hapalt = copy.copy(self.vrnt_hapalt),
            vrnt_hapref = copy.copy(self.vrnt_hapref),
            vrnt_mask = copy.copy(self.vrnt_mask)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len = copy.copy(self.vrnt_chrgrp_len)

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Deep copy metadata.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A deep copy of the DenseTaxaVariantMatrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            vrnt_chrgrp = copy.deepcopy(self.vrnt_chrgrp, memo),
            vrnt_phypos = copy.deepcopy(self.vrnt_phypos, memo),
            vrnt_name = copy.deepcopy(self.vrnt_name, memo),
            vrnt_genpos = copy.deepcopy(self.vrnt_genpos, memo),
            vrnt_xoprob = copy.deepcopy(self.vrnt_xoprob, memo),
            vrnt_hapgrp = copy.deepcopy(self.vrnt_hapgrp, memo),
            vrnt_hapalt = copy.deepcopy(self.vrnt_hapalt, memo),
            vrnt_hapref = copy.deepcopy(self.vrnt_hapref, memo),
            vrnt_mask = copy.deepcopy(self.vrnt_mask, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name, memo)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix, memo)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix, memo)
        out.vrnt_chrgrp_len = copy.deepcopy(self.vrnt_chrgrp_len, memo)

        return out

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseTaxaMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set raw underlying numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim_gteq(value, "mat", 2)
        self._mat = value

    ############### Taxa Metadata Properites ###############
    @DenseTaxaMatrix.taxa_axis.getter
    def taxa_axis(self) -> int:
        """Get taxa axis number"""
        return 0

    ############# Variant Metadata Properites ##############
    @DenseVariantMatrix.vrnt_axis.getter
    def vrnt_axis(self) -> int:
        """Get variant axis"""
        return 1

    ############################## Object Methods ##############################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapalt : numpy.ndarray
            Variant alternative haplotype labels to adjoin to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapref : numpy.ndarray
            Variant reference haplotype labels to adjoin to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A copy of DenseTaxaVariantMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.adjoin_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
            out = self.adjoin_vrnt(
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A copy of DenseTaxaVariantMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseTaxaVariantMatrix is allocated and filled.

        """
        out = super(DenseTaxaVariantMatrix, self).adjoin_taxa(
            values = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def adjoin_vrnt(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Add additional elements to the end of the Matrix along the variant axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to adjoin to the Matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to adjoin to the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to adjoin to the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to adjoin to the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to adjoin to the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to adjoin to the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to adjoin to the Matrix.
        vrnt_hapalt : numpy.ndarray
            Variant haplotype sequence.
        vrnt_hapref : numpy.ndarray
            Variant haplotype reference sequence.
        vrnt_mask : numpy.ndarray
            Variant mask to adjoin to the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A copy of DenseTaxaVariantMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, self).adjoin_vrnt(
            values = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len

        return out

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Delete sub-arrays along an axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to delete the subarray defined by obj.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A DenseTaxaVariantMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.delete_taxa(
                obj = obj,
                **kwargs
            )
        elif axis == self.vrnt_axis:
            out = self.delete_vrnt(
                obj = obj,
                **kwargs
            )
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Delete sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A DenseTaxaVariantMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, self).delete_taxa(
            obj = obj,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def delete_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Delete sub-arrays along the variant axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A DenseTaxaVariantMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, self).delete_vrnt(
            obj = obj,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len

        return out

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapalt : numpy.ndarray
            Variant alternative haplotype labels to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapref : numpy.ndarray
            Variant reference haplotype labels to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A DenseTaxaVariantMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.insert_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
            out = self.insert_vrnt(
                obj = obj,
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    def insert_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Insert values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A DenseTaxaVariantMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        # create output
        out = super(DenseTaxaVariantMatrix, self).insert_taxa(
            obj = obj,
            values = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def insert_vrnt(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Insert values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to insert into the Matrix.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to insert into the Matrix.
        vrnt_name : numpy.ndarray
            Variant names to insert into the Matrix.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to insert into the Matrix.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to insert into the Matrix.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to insert into the Matrix.
        vrnt_hapalt : numpy.ndarray
            Variant alternative haplotype labels to insert into the Matrix.
        vrnt_hapref : numpy.ndarray
            Variant reference haplotype labels to insert into the Matrix.
        vrnt_mask : numpy.ndarray
            Variant mask to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            A DenseTaxaVariantMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        # create output
        out = super(DenseTaxaVariantMatrix, self).insert_vrnt(
            obj = obj,
            values = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len

        return out

    def select(
            self, 
            indices: ArrayLike, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Select certain values from the matrix.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        axis : int
            The axis along which values are selected.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            The output DenseTaxaVariantMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.select_vrnt(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_taxa(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Select certain values from the Matrix along the taxa axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            The output DenseTaxaVariantMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, self).select_taxa(
            indices = indices,
            vrnt_chrgrp = self._vrnt_chrgrp,
            vrnt_phypos = self._vrnt_phypos,
            vrnt_name = self._vrnt_name,
            vrnt_genpos = self._vrnt_genpos,
            vrnt_xoprob = self._vrnt_xoprob,
            vrnt_hapgrp = self._vrnt_hapgrp,
            vrnt_hapalt = self._vrnt_hapalt,
            vrnt_hapref = self._vrnt_hapref,
            vrnt_mask = self._vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def select_vrnt(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Select certain values from the Matrix along the variant axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaVariantMatrix
            The output DenseTaxaVariantMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaVariantMatrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, self).select_vrnt(
            indices = indices,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len

        return out

    @classmethod
    def concat(
            cls, 
            mats: Sequence, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Concatenate matrices together along an axis.

        Parameters
        ----------
        mats : Sequence of matrices
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : DenseTaxaVariantMatrix
            The concatenated DenseTaxaVariantMatrix. Note that concat does not occur in-place:
            a new DenseTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].taxa_axis:
            out = cls.concat_taxa(mats, **kwargs)
        elif axis == mats[0].vrnt_axis:
            out = cls.concat_vrnt(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_taxa(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Concatenate list of Matrix together along the taxa axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : DenseTaxaVariantMatrix
            The concatenated DenseTaxaVariantMatrix. Note that concat does not occur in-place:
            a new DenseTaxaVariantMatrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, cls).concat_taxa(
            mats = mats,
            vrnt_chrgrp = mats[0].vrnt_chrgrp,
            vrnt_phypos = mats[0].vrnt_phypos,
            vrnt_name = mats[0].vrnt_name,
            vrnt_genpos = mats[0].vrnt_genpos,
            vrnt_xoprob = mats[0].vrnt_xoprob,
            vrnt_hapgrp = mats[0].vrnt_hapgrp,
            vrnt_hapalt = mats[0].vrnt_hapalt,
            vrnt_hapref = mats[0].vrnt_hapref,
            vrnt_mask = mats[0].vrnt_mask,
            **kwargs
        )

        # copy metadata from source
        out.vrnt_chrgrp_name = mats[0].vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = mats[0].vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = mats[0].vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = mats[0].vrnt_chrgrp_len

        return out

    @classmethod
    def concat_vrnt(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseTaxaVariantMatrix':
        """
        Concatenate list of Matrix together along the variant axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        out = super(DenseTaxaVariantMatrix, cls).concat_vrnt(
            mats = mats,
            taxa = mats[0].taxa,
            taxa_grp = mats[0].taxa_grp,
            **kwargs
        )

        # copy metadata from source
        out.taxa_grp_name = mats[0].taxa_grp_name
        out.taxa_grp_stix = mats[0].taxa_grp_stix
        out.taxa_grp_spix = mats[0].taxa_grp_spix
        out.taxa_grp_len = mats[0].taxa_grp_len

        return out

    ######### Matrix element in-place-manipulation #########
    def append(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseHaplotypeMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapalt : numpy.ndarray
            Variant alternative haplotype labels to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapref : numpy.ndarray
            Variant reference haplotype labels to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to append to the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.append_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
            self.append_vrnt(
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    def remove(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along an axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.taxa_axis:
            self.remove_taxa(obj = obj, **kwargs)
        elif axis == self.vrnt_axis:
            self.remove_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : array_like
            Values to incorporate into the matrix.
        axis : int
            The axis along which values are incorporated.
        taxa : numpy.ndarray
            Taxa names to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        vrnt_chrgrp : numpy.ndarray
            Variant chromosome groups to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_chrgrp field, providing this argument overwrites the field.
        vrnt_phypos : numpy.ndarray
            Variant chromosome physical positions to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_phypos field, providing this argument overwrites the field.
        vrnt_name : numpy.ndarray
            Variant names to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_name field, providing this argument overwrites the field.
        vrnt_genpos : numpy.ndarray
            Variant chromosome genetic positions to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_genpos field, providing this argument overwrites the field.
        vrnt_xoprob : numpy.ndarray
            Sequential variant crossover probabilities to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_xoprob field, providing this argument overwrites the field.
        vrnt_hapgrp : numpy.ndarray
            Variant haplotype labels to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapalt : numpy.ndarray
            Variant alternative haplotype labels to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_hapref : numpy.ndarray
            Variant reference haplotype labels to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_hapgrp field, providing this argument overwrites the field.
        vrnt_mask : numpy.ndarray
            Variant mask to incorporate into the Matrix.
            If values is a DenseTaxaVariantMatrix that has a non-None
            vrnt_mask field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.taxa_axis:
            self.incorp_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis == self.vrnt_axis:
            self.incorp_vrnt(
                obj = obj,
                values = values,
                vrnt_chrgrp = vrnt_chrgrp,
                vrnt_phypos = vrnt_phypos,
                vrnt_name = vrnt_name,
                vrnt_genpos = vrnt_genpos,
                vrnt_xoprob = vrnt_xoprob,
                vrnt_hapgrp = vrnt_hapgrp,
                vrnt_hapalt = vrnt_hapalt,
                vrnt_hapref = vrnt_hapref,
                vrnt_mask = vrnt_mask,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    ################### Sorting Methods ####################
    def lexsort(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            axis: int = -1, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            The axis of the Matrix over which to sort values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
        """
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        indices = None                          # declare variable

        # dispatch to correct function
        if axis == self.taxa_axis:
            indices = self.lexsort_taxa(keys = keys, **kwargs)
        elif axis == self.vrnt_axis:
            indices = self.lexsort_vrnt(keys = keys, **kwargs)
        else:
            raise ValueError("cannot lexsort along axis {0}".format(axis))

        return indices

    def reorder(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Reorder the VariantMatrix.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            The axis over which to reorder values.
        kwargs : dict
            Additional keyword arguments.
        """
        axis = get_axis(axis, self.mat_ndim)                   # transform axis number to an index

        if axis == self.taxa_axis:
            self.reorder_taxa(indices = indices, **kwargs)
        elif axis == self.vrnt_axis:
            self.reorder_vrnt(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def sort(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Reset metadata for corresponding axis: name, stix, spix, len.
        Sort the VariantMatrix using a tuple of keys.

        Parameters
        ----------
        keys : tuple, None
            A tuple of columns to be sorted. The last column is the primary
            sort key. If None, sort using vrnt_chrgrp as primary key, and
            vrnt_phypos as secondary key.
        axis : int
            The axis over which to sort values.
        kwargs : dict
            Additional keyword arguments.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.sort_taxa(keys = keys, **kwargs)
        elif axis == self.vrnt_axis:
            self.sort_vrnt(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    ################### Grouping Methods ###################
    def group(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Sort the DenseTaxaVariantMatrix along an axis, then populate grouping 
        indices.

        Parameters
        ----------
        axis : int
            The axis along which values are grouped.
        kwargs : dict
            Additional keyword arguments.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.group_taxa(**kwargs)
        elif axis == self.vrnt_axis:
            self.group_vrnt(**kwargs)
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

    def ungroup(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Ungroup the DenseTaxaVariantMatrix along an axis by removing grouping 
        metadata.

        Parameters
        ----------
        axis : int
            The axis along which values should be ungrouped.
        kwargs : dict
            Additional keyword arguments.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.ungroup_taxa(**kwargs)
        elif axis == self.vrnt_axis:
            self.ungroup_vrnt(**kwargs)
        else:
            raise ValueError("cannot ungroup along axis {0}".format(axis))

    def is_grouped(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the Matrix has been sorted and grouped.

        Parameters
        ----------
        axis : int
            Axis along which to determine if is grouped.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        grouped = False                         # default output

        if axis == self.taxa_axis:
            grouped = self.is_grouped_taxa(**kwargs)
        elif axis == self.vrnt_axis:
            grouped = self.is_grouped_vrnt(**kwargs)
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write ``DenseTaxaVariantMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseTaxaVariantMatrix`` data is stored.
            If None, ``DenseTaxaVariantMatrix`` is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite values in an HDF5 file if a field already exists.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r+``) mode
        if isinstance(filename, (str,Path)):
            h5file = h5py.File(filename, "a")

        # elif we have an h5py.File, make sure mode is writable, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_writable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        #### write data to HDF5 file and (optionally) close ####

        # data dictionary
        data = {
            "mat"               : self.mat,
            "taxa"              : self.taxa,
            "taxa_grp"          : self.taxa_grp,
            "vrnt_chrgrp"       : self.vrnt_chrgrp,
            "vrnt_phypos"       : self.vrnt_phypos,
            "vrnt_name"         : self.vrnt_name,
            "vrnt_genpos"       : self.vrnt_genpos,
            "vrnt_xoprob"       : self.vrnt_xoprob,
            "vrnt_hapgrp"       : self.vrnt_hapgrp,
            "vrnt_hapalt"       : self.vrnt_hapalt,
            "vrnt_hapref"       : self.vrnt_hapref,
            "vrnt_mask"         : self.vrnt_mask,
            # metadata
            "taxa_grp_name"     : self.taxa_grp_name,
            "taxa_grp_stix"     : self.taxa_grp_stix,
            "taxa_grp_spix"     : self.taxa_grp_spix,
            "taxa_grp_len"      : self.taxa_grp_len,
            "vrnt_chrgrp_name"  : self.vrnt_chrgrp_name,
            "vrnt_chrgrp_stix"  : self.vrnt_chrgrp_stix,
            "vrnt_chrgrp_spix"  : self.vrnt_chrgrp_spix,
            "vrnt_chrgrp_len"   : self.vrnt_chrgrp_len,
        }

        # save data
        h5py_File_write_dict(h5file, groupname, data, overwrite)

        # close the file, only if the provided filename was a string or Path and not a h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseTaxaVariantMatrix':
        """
        Read ``DenseTaxaVariantMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseTaxaVariantMatrix`` data is stored.
            If None, ``DenseTaxaVariantMatrix`` is read from base HDF5 group.

        Returns
        -------
        gmat : DenseTaxaVariantMatrix
            A genotype matrix read from file.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in read (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an ``h5py.File``, make sure mode is in at least ``r`` mode, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_readable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # FIXME: errors if groupname == "" or "/"
            # if the group does not exist in the file, close and raise error
            check_h5py_File_has_group(h5file, groupname)

            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = ["mat"]

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)
        
        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
            "mat"               : None,
            "taxa"              : None,
            "taxa_grp"          : None,
            "vrnt_chrgrp"       : None,
            "vrnt_phypos"       : None,
            "vrnt_name"         : None,
            "vrnt_genpos"       : None,
            "vrnt_xoprob"       : None,
            "vrnt_hapgrp"       : None,
            "vrnt_hapalt"       : None,
            "vrnt_hapref"       : None,
            "vrnt_mask"         : None,
            # metadata
            "taxa_grp_name"     : None,
            "taxa_grp_stix"     : None,
            "taxa_grp_spix"     : None,
            "taxa_grp_len"      : None,
            "vrnt_chrgrp_name"  : None,
            "vrnt_chrgrp_stix"  : None,
            "vrnt_chrgrp_spix"  : None,
            "vrnt_chrgrp_len"   : None,
        }

        ##################################
        ### read mandatory data fields ###

        # read mat array (ndarray dtype = any)
        data["mat"] = h5py_File_read_ndarray(h5file, groupname + "mat")
        
        #################################
        ### read optional data fields ###

        # read taxa array (ndarray dtype = unicode / object)
        if groupname + "taxa" in h5file:
            data["taxa"] = h5py_File_read_ndarray_utf8(h5file, groupname + "taxa")

        # read taxa_grp array (ndarray dtype = any)
        if groupname + "taxa_grp" in h5file:
            data["taxa_grp"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp")
        
        # read vrnt_chrgrp array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp" in h5file:
            data["vrnt_chrgrp"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp")

        # read vrnt_phypos array (ndarray dtype = any)
        if groupname + "vrnt_phypos" in h5file:
            data["vrnt_phypos"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_phypos")

        # read vrnt_name array (ndarray dtype = unicode / object)
        if groupname + "vrnt_name" in h5file:
            data["vrnt_name"] = h5py_File_read_ndarray_utf8(h5file, groupname + "vrnt_name")

        # read vrnt_genpos array (ndarray dtype = any)
        if groupname + "vrnt_genpos" in h5file:
            data["vrnt_genpos"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_genpos")

        # read vrnt_xoprob array (ndarray dtype = any)
        if groupname + "vrnt_xoprob" in h5file:
            data["vrnt_xoprob"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_xoprob")

        # read vrnt_hapgrp array (ndarray dtype = any)
        if groupname + "vrnt_hapgrp" in h5file:
            data["vrnt_hapgrp"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_hapgrp")

        # read vrnt_hapalt array (ndarray dtype = unicode / object)
        if groupname + "vrnt_hapalt" in h5file:
            data["vrnt_hapalt"] = h5py_File_read_ndarray_utf8(h5file, groupname + "vrnt_hapalt")

        # read vrnt_hapref array (ndarray dtype = unicode / object)
        if groupname + "vrnt_hapref" in h5file:
            data["vrnt_hapref"] = h5py_File_read_ndarray_utf8(h5file, groupname + "vrnt_hapref")

        # read vrnt_mask array (ndarray dtype = any)
        if groupname + "vrnt_mask" in h5file:
            data["vrnt_mask"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_mask")

        #####################################
        ### read optional metadata fields ###

        # read taxa_grp_name array (ndarray dtype = any)
        if groupname + "taxa_grp_name" in h5file:
            data["taxa_grp_name"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_name")

        # read taxa_grp_stix array (ndarray dtype = any)
        if groupname + "taxa_grp_stix" in h5file:
            data["taxa_grp_stix"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_stix")

        # read taxa_grp_spix array (ndarray dtype = any)
        if groupname + "taxa_grp_spix" in h5file:
            data["taxa_grp_spix"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_spix")

        # read taxa_grp_len array (ndarray dtype = any)
        if groupname + "taxa_grp_len" in h5file:
            data["taxa_grp_len"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp_len")

        # read vrnt_chrgrp_name array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_name" in h5file:
            data["vrnt_chrgrp_name"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_name")

        # read vrnt_chrgrp_stix array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_stix" in h5file:
            data["vrnt_chrgrp_stix"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_stix")

        # read vrnt_chrgrp_spix array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_spix" in h5file:
            data["vrnt_chrgrp_spix"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_spix")

        # read vrnt_chrgrp_len array (ndarray dtype = any)
        if groupname + "vrnt_chrgrp_len" in h5file:
            data["vrnt_chrgrp_len"] = h5py_File_read_ndarray(h5file, groupname + "vrnt_chrgrp_len")

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ################### Object creation ####################
        
        # create object from read data
        out = cls(
            mat         = data["mat"],
            taxa        = data["taxa"],
            taxa_grp    = data["taxa_grp"],
            vrnt_chrgrp = data["vrnt_chrgrp"],
            vrnt_phypos = data["vrnt_phypos"],
            vrnt_name   = data["vrnt_name"],
            vrnt_genpos = data["vrnt_genpos"],
            vrnt_xoprob = data["vrnt_xoprob"],
            vrnt_hapgrp = data["vrnt_hapgrp"],
            vrnt_hapalt = data["vrnt_hapalt"],
            vrnt_hapref = data["vrnt_hapref"],
            vrnt_mask   = data["vrnt_mask"], 
        )

        # copy metadata
        out.taxa_grp_name    = data["taxa_grp_name"]
        out.taxa_grp_stix    = data["taxa_grp_stix"]
        out.taxa_grp_spix    = data["taxa_grp_spix"]
        out.taxa_grp_len     = data["taxa_grp_len"]
        out.vrnt_chrgrp_name = data["vrnt_chrgrp_name"]
        out.vrnt_chrgrp_stix = data["vrnt_chrgrp_stix"]
        out.vrnt_chrgrp_spix = data["vrnt_chrgrp_spix"]
        out.vrnt_chrgrp_len  = data["vrnt_chrgrp_len"]

        return out



################################## Utilities ###################################
def check_is_DenseTaxaVariantMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseTaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseTaxaVariantMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseTaxaVariantMatrix.__name__,type(v).__name__))
