"""
Module implementing a dense matrix with phase, taxa, and variant metadata
and associated error checking routines.
"""

import numpy
from typing import Any, Sequence, Union
from typing import Optional
from numpy.typing import ArrayLike

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_ndim_gteq
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.core.mat.DensePhasedMatrix import DensePhasedMatrix
from pybrops.core.mat.DenseTaxaVariantMatrix import DenseTaxaVariantMatrix

class DensePhasedTaxaVariantMatrix(DenseTaxaVariantMatrix,DensePhasedMatrix,PhasedTaxaVariantMatrix):
    """
    A concrete class implementing dense matrices with phase, taxa, and variant
    metadata.

    The purpose of this concrete class is to merge the following classes and
    interfaces:

        1. DenseTaxaVariantMatrix
        2. DensePhasedMatrix
        3. PhasedTaxaVariantMatrix
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
        Parameters
        ----------
        mat : numpy.ndarray
            An int8 haplotype matrix. Must be {0,1,2} format.
        taxa : numpy.ndarray, None
        taxa_grp : numpy.ndarray, None
        vrnt_chrgrp : numpy.ndarray, None
        vrnt_phypos : numpy.ndarray, None
        vrnt_name : numpy.ndarray, None
        vrnt_genpos : numpy.ndarray, None
        vrnt_xoprob : numpy.ndarray, None
        vrnt_hapgrp : numpy.ndarray, None
        vrnt_hapalt : numpy.ndarray, None
        vrnt_hapref : numpy.ndarray, None
        vrnt_mask : numpy.ndarray, None
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DensePhasedTaxaVariantMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
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

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ##################### Matrix Data ######################
    @DenseTaxaVariantMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        check_is_ndarray(value, "mat")
        check_ndarray_ndim_gteq(value, "mat", 3)
        self._mat = value

    ############## Phase Metadata Properites ###############
    @DensePhasedMatrix.phase_axis.getter
    def phase_axis(self) -> int:
        """Get phase axis number"""
        return 0

    ############### Taxa Metadata Properites ###############
    @DenseTaxaVariantMatrix.taxa_axis.getter
    def taxa_axis(self) -> int:
        """Get taxa axis number"""
        return 1

    ############# Variant Metadata Properites ##############
    @DenseTaxaVariantMatrix.vrnt_axis.getter
    def vrnt_axis(self) -> int:
        """Get variant axis"""
        return 2

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
        ) -> 'DensePhasedTaxaVariantMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedTaxaVariantMatrix
            A copy of DensePhasedTaxaVariantMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.adjoin_phase(
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
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

    def adjoin_phase(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
        """
        Adjoin values along the phase axis.

        Parameters
        ----------
        values : Matrix or numpy.ndarray
            Values to adjoin along the phase axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedTaxaVariantMatrix
            A copy of DensePhasedTaxaVariantMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        out = super(DensePhasedTaxaVariantMatrix, self).adjoin_phase(
            values = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
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
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
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
        out : DensePhasedTaxaVariantMatrix
            A DensePhasedTaxaVariantMatrix with deleted elements. Note that concat does not occur
            in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.delete_phase(obj = obj, **kwargs)
        elif axis == self.taxa_axis:
            out = self.delete_taxa(obj = obj, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.delete_vrnt(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
        """
        Delete sub-arrays along the phase axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedTaxaVariantMatrix
            A DensePhasedTaxaVariantMatrix with deleted elements. Note that concat does not occur
            in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        out = super(DensePhasedTaxaVariantMatrix, self).delete_phase(
            obj = obj,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
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
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

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
        ) -> 'DensePhasedTaxaVariantMatrix':
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values: Matrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedTaxaVariantMatrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.insert_phase(
                obj = obj,
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
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

    def insert_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
        """
        Insert values along the phase axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedTaxaVariantMatrix
            A DensePhasedTaxaVariantMatrix with values inserted. Note that insert does not occur
            in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        out = super(DensePhasedTaxaVariantMatrix, self).insert_phase(
            obj = obj,
            values = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
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
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    def select(
            self, 
            indices: ArrayLike, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
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
        out : DensePhasedTaxaVariantMatrix
            The output DensePhasedTaxaVariantMatrix with values selected. Note that select does not
            occur in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.select_phase(indices = indices, **kwargs)
        elif axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis == self.vrnt_axis:
            out = self.select_vrnt(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_phase(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
        """
        Select certain values from the Matrix along the phase axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedTaxaVariantMatrix
            The output DensePhasedTaxaVariantMatrix with values selected. Note that select does not
            occur in-place: a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        out = super(DensePhasedTaxaVariantMatrix, self).select_phase(
            indices = indices,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
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
        out.taxa_grp_name = self._taxa_grp_name
        out.taxa_grp_stix = self._taxa_grp_stix
        out.taxa_grp_spix = self._taxa_grp_spix
        out.taxa_grp_len = self._taxa_grp_len
        out.vrnt_chrgrp_name = self._vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = self._vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = self._vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = self._vrnt_chrgrp_len

        return out

    @classmethod
    def concat(
            cls, 
            mats: Sequence, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
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
        out : DensePhasedTaxaVariantMatrix
            The concatenated DensePhasedTaxaVariantMatrix. Note that concat does not occur in-place:
            a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].phase_axis:
            out = cls.concat_phase(mats, **kwargs)
        elif axis == mats[0].taxa_axis:
            out = cls.concat_taxa(mats, **kwargs)
        elif axis == mats[0].vrnt_axis:
            out = cls.concat_vrnt(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_phase(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DensePhasedTaxaVariantMatrix':
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
        out : DensePhasedTaxaVariantMatrix
            The concatenated DensePhasedTaxaVariantMatrix. Note that concat does not occur in-place:
            a new DensePhasedTaxaVariantMatrix is allocated and filled.
        """
        out = super(DensePhasedTaxaVariantMatrix, cls).concat_phase(
            mats = mats,
            taxa = mats[0].taxa,
            taxa_grp = mats[0].taxa_grp,
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
        out.taxa_grp_name = mats[0].taxa_grp_name
        out.taxa_grp_stix = mats[0].taxa_grp_stix
        out.taxa_grp_spix = mats[0].taxa_grp_spix
        out.taxa_grp_len = mats[0].taxa_grp_len
        out.vrnt_chrgrp_name = mats[0].vrnt_chrgrp_name
        out.vrnt_chrgrp_stix = mats[0].vrnt_chrgrp_stix
        out.vrnt_chrgrp_spix = mats[0].vrnt_chrgrp_spix
        out.vrnt_chrgrp_len = mats[0].vrnt_chrgrp_len

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
        values: Matrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.phase_axis:
            self.append_phase(
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
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

        if axis == self.phase_axis:
            self.remove_phase(obj = obj, **kwargs)
        elif axis == self.taxa_axis:
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
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.phase_axis:
            self.incorp_phase(
                obj = obj,
                values = values,
                **kwargs
            )
        elif axis == self.taxa_axis:
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

    ################### Grouping Methods ###################
    def is_grouped(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the Matrix has been sorted and grouped.

        Returns
        -------
        grouped : bool
            True or False indicating whether the GeneticMap has been sorted and
            grouped.
        """
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        grouped = False                         # default output

        if axis == self.phase_axis:
            grouped = False     # always results to false, even though axis is ungroupable
        elif axis == self.taxa_axis:
            grouped = self.is_grouped_taxa(**kwargs)
        elif axis == self.vrnt_axis:
            grouped = self.is_grouped_vrnt(**kwargs)
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped



################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedTaxaVariantMatrix(v: object) -> bool:
    """
    Determine whether an object is a DensePhasedTaxaVariantMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DensePhasedTaxaVariantMatrix object instance.
    """
    return isinstance(v, DensePhasedTaxaVariantMatrix)

def check_is_DensePhasedTaxaVariantMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DensePhasedTaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DensePhasedTaxaVariantMatrix):
        raise TypeError("'{0}' must be a DensePhasedTaxaVariantMatrix".format(vname))
