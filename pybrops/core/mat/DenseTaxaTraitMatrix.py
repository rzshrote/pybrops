"""
Module defining implementing dense matrices with taxa and trait metadata and
associated error checking routines.
"""

__all__ = [
    "DenseTaxaTraitMatrix",
    "check_is_DenseTaxaTraitMatrix"
]

import numpy
from numpy.typing import ArrayLike
import copy
from typing import Optional, Sequence, Union

from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.mat.DenseTaxaMatrix import DenseTaxaMatrix
from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.core.mat.TaxaTraitMatrix import TaxaTraitMatrix

class DenseTaxaTraitMatrix(DenseTaxaMatrix,DenseTraitMatrix,TaxaTraitMatrix):
    """
    A concrete class for dense matrices with taxa and trait metadata.

    The purpose of this concrete class is to merge the following implementations
    and interfaces:

        1. DenseTaxaMatrix
        2. DenseTraitMatrix
        3. TaxaTraitMatrix
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
        Constructor for the concrete class DenseTaxaTraitMatrix.

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
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
            Additional keyword arguments.
        """
        # cannot use super constructor since DenseTaxaTraitMatrix and DenseVariantMatrix
        # constructors have overlapping arguments

        # set data
        self.mat = mat
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.trait = trait

        # set taxa metadata
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseTaxaTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait)
        )
        # copy taxa metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseTaxaTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        return out

    ############################ Object Properties #############################

    ############### Taxa Metadata Properites ###############
    @DenseTaxaMatrix.taxa_axis.getter
    def taxa_axis(self) -> int:
        """Get taxa axis number"""
        return 0

    ############# Variant Metadata Properites ##############
    @DenseTraitMatrix.trait_axis.getter
    def trait_axis(self):
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
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DenseTaxaTraitMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Traits to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaTraitMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
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
        elif axis == self.trait_axis:
            out = self.adjoin_trait(
                values = values,
                trait = trait,
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
        ) -> 'DenseTaxaTraitMatrix':
        """
        Add additional elements to the end of the Matrix along the taxa axis.

        Parameters
        ----------
        values : DenseTaxaTraitMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Traits to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaTraitMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, self).adjoin_taxa(
            values = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self._trait,
            **kwargs
        )

        return out

    def adjoin_trait(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
        """
        Add additional elements to the end of the Matrix along the trait axis.

        Parameters
        ----------
        values : DenseTaxaTraitMatrix, numpy.ndarray
            Values to be adjoined to the Matrix.
        trait : numpy.ndarray
            Trait names to adjoin to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaTraitMatrix
            A copy of mat with values adjoined to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, self).adjoin_trait(
            values = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            trait = trait,
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
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            A DenseTaxaTraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.delete_taxa(
                obj = obj,
                **kwargs
            )
        elif axis == self.trait_axis:
            out = self.delete_trait(
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
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            A DenseTaxaTraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, self).delete_taxa(
            obj = obj,
            trait = self._trait,
            **kwargs
        )

        return out

    def delete_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            A DenseTaxaTraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, self).delete_trait(
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
            axis: int = -1, taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : DenseTaxaTraitMatrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        taxa : numpy.ndarray
            Taxa names to insert into the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to insert into the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaTraitMatrix
            A DenseTaxaTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaTraitMatrix is allocated and filled.
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
        elif axis == self.trait_axis:
            out = self.insert_trait(
                obj = obj,
                values = values,
                trait = trait,
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
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            A DenseTaxaTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        # create output
        out = super(DenseTaxaTraitMatrix, self).insert_taxa(
            obj = obj,
            values = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self._trait,
            **kwargs
        )

        return out

    def insert_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
        """
        Insert values along the variant axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : array_like
            Values to insert into the matrix.
        trait : numpy.ndarray
            Trait names to insert into the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaTraitMatrix
            A DenseTaxaTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        # create output
        out = super(DenseTaxaTraitMatrix, self).insert_trait(
            obj = obj,
            values = values,
            taxa = self._taxa,
            taxa_grp = self._taxa_grp,
            trait = trait,
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
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            The output DenseTaxaTraitMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis == self.trait_axis:
            out = self.select_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_taxa(
            self, 
            indices, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            The output DenseTaxaTraitMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, self).select_taxa(
            indices = indices,
            trait = self._trait,
            **kwargs
        )

        return out

    def select_trait(
            self, 
            indices, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            The output DenseTaxaTraitMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaTraitMatrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, self).select_trait(
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
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].taxa_axis:
            out = cls.concat_taxa(mats, **kwargs)
        elif axis == mats[0].trait_axis:
            out = cls.concat_trait(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_taxa(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            The concatenated DenseTaxaTraitMatrix. Note that concat does not occur in-place:
            a new DenseTaxaTraitMatrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, cls).concat_taxa(
            mats = mats,
            trait = mats[0].trait,
            **kwargs
        )

        return out

    @classmethod
    def concat_trait(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseTaxaTraitMatrix':
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
        out : DenseTaxaTraitMatrix
            The concatenated DenseTaxaTraitMatrix. Note that concat does not occur in-place:
            a new DenseTaxaTraitMatrix is allocated and filled.
        """
        out = super(DenseTaxaTraitMatrix, cls).concat_trait(
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
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseTaxaTraitMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to append to the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
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
        elif axis == self.trait_axis:
            self.append_trait(
                values = values,
                trait = trait,
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
        elif axis == self.trait_axis:
            self.remove_trait(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
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
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to incorporate into the Matrix.
            If values is a DenseTaxaTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
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
        elif axis == self.trait_axis:
            self.incorp_trait(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    ################### Sorting Methods ####################
    def lexsort(
            self, 
            keys: Union[tuple,numpy.ndarray,None], 
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

        Returns
        -------
        indices : numpy.ndarray
            Array of indices that sort the keys.
        """
        axis = get_axis(axis, self.mat_ndim)    # transform axis number to an index
        indices = None                          # declare variable

        # dispatch to correct function
        if axis == self.taxa_axis:
            self.lexsort_taxa(keys = keys, **kwargs)
        elif axis == self.trait_axis:
            self.lexsort_trait(keys = keys, **kwargs)
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
        """
        axis = get_axis(axis, self.mat_ndim)                   # transform axis number to an index

        if axis == self.taxa_axis:
            self.reorder_taxa(indices = indices, **kwargs)
        elif axis == self.trait_axis:
            self.reorder_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def sort(
            self, 
            keys: Union[tuple,numpy.ndarray,None], 
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
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.sort_taxa(keys = keys, **kwargs)
        elif axis == self.trait_axis:
            self.sort_trait(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    ################### Grouping Methods ###################
    def group(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        Calculate chromosome grouping indices (group by vrnt_chrgrp).
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.taxa_axis:
            self.group_taxa(**kwargs)
        elif axis == self.trait_axis:
            raise ValueError("cannot group along axis {0} (trait axis)".format(axis))
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

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

        if axis == self.taxa_axis:
            grouped = self.is_grouped_taxa(**kwargs)
        elif axis == self.trait_axis:
            raise ValueError("cannot test for grouping along axis {0} (trait axis)".format(axis))
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped



################################## Utilities ###################################
def check_is_DenseTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseTaxaTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseTaxaTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseTaxaTraitMatrix.__name__,type(v).__name__))
