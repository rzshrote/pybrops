"""
Module defining implementing dense matrices with taxa metadata and associated
error checking routines.
"""

__all__ = [
    "DenseTaxaMatrix",
    "check_is_DenseTaxaMatrix",
]

import copy
from pathlib import Path
import numpy
import h5py
from typing import Optional
from typing import Sequence
from typing import Union

from numpy.typing import ArrayLike

from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_python import check_is_array_like
from pybrops.core.error.error_attr_python import check_is_iterable
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_integer
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_generic_python import generic_check_isinstance
from pybrops.core.error.error_value_numpy import check_ndarray_ndim_gteq
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.util.array import get_axis
from pybrops.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.util.h5py import h5py_File_read_ndarray
from pybrops.core.util.h5py import h5py_File_read_ndarray_utf8
from pybrops.core.util.h5py import h5py_File_write_dict

class DenseTaxaMatrix(
        DenseMutableMatrix,
        TaxaMatrix,
    ):
    """
    A concrete class for dense matrices with taxa metadata.

    The purpose of this concrete class is to implement base functionality for:
        1) Dense matrix taxa metadata.
        2) Dense matrix taxa routines.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseTaxaMatrix.

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
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseTaxaMatrix, self).__init__(
            mat = mat,
            **kwargs
        )
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseTaxaMatrix':
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
            taxa_grp = copy.copy(self.taxa_grp)
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
        ) -> 'DenseTaxaMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        return out

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    @DenseMutableMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        """Set pointer to raw numpy.ndarray object."""
        check_is_ndarray(value, "mat")
        check_ndarray_ndim_gteq(value, "mat", 1)
        self._mat = value

    ################# Taxa Data Properites #################
    @property
    def taxa(self) -> Union[numpy.ndarray,None]:
        """Taxa label array"""
        return self._taxa
    @taxa.setter
    def taxa(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa label array"""
        if value is not None:
            check_is_ndarray(value, "taxa")
            check_ndarray_dtype_is_object(value, "taxa")
            check_ndarray_ndim(value, "taxa", 1)
            check_ndarray_axis_len(value, "taxa", 0, self.ntaxa)
        self._taxa = value

    @property
    def taxa_grp(self) -> Union[numpy.ndarray,None]:
        """Taxa group label."""
        return self._taxa_grp
    @taxa_grp.setter
    def taxa_grp(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group label array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp")
            check_ndarray_dtype_is_integer(value, "taxa_grp")
            check_ndarray_ndim(value, "taxa_grp", 1)
            check_ndarray_axis_len(value, "taxa_grp", 0, self.ntaxa)
        self._taxa_grp = value

    ############### Taxa Metadata Properites ###############
    @property
    def ntaxa(self) -> int:
        """Number of taxa"""
        return self._mat.shape[self.taxa_axis]
    
    @property
    def taxa_axis(self) -> int:
        """Axis along which taxa are stored"""
        return 0
    
    @property
    def taxa_grp_name(self) -> Union[numpy.ndarray,None]:
        """Taxa group name."""
        return self._taxa_grp_name
    @taxa_grp_name.setter
    def taxa_grp_name(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group name array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_name")
            check_ndarray_dtype_is_integer(value, "taxa_grp_name")
            check_ndarray_ndim(value, "taxa_grp_name", 1)
        self._taxa_grp_name = value

    @property
    def taxa_grp_stix(self) -> Union[numpy.ndarray,None]:
        """Taxa group start index."""
        return self._taxa_grp_stix
    @taxa_grp_stix.setter
    def taxa_grp_stix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group start indices array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_stix")
            check_ndarray_dtype_is_integer(value, "taxa_grp_stix")
            check_ndarray_ndim(value, "taxa_grp_stix", 1)
        self._taxa_grp_stix = value

    @property
    def taxa_grp_spix(self) -> Union[numpy.ndarray,None]:
        """Taxa group stop index."""
        return self._taxa_grp_spix
    @taxa_grp_spix.setter
    def taxa_grp_spix(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group stop indices array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_spix")
            check_ndarray_dtype_is_integer(value, "taxa_grp_spix")
            check_ndarray_ndim(value, "taxa_grp_spix", 1)
        self._taxa_grp_spix = value

    @property
    def taxa_grp_len(self) -> Union[numpy.ndarray,None]:
        """Taxa group length."""
        return self._taxa_grp_len
    @taxa_grp_len.setter
    def taxa_grp_len(self, value: Union[numpy.ndarray,None]) -> None:
        """Set taxa group length array"""
        if value is not None:
            check_is_ndarray(value, "taxa_grp_len")
            check_ndarray_dtype_is_integer(value, "taxa_grp_len")
            check_ndarray_ndim(value, "taxa_grp_len", 1)
        self._taxa_grp_len = value

    ############################## Object Methods ##############################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values to be appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaMatrix
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
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
        """
        Add additional elements to the end of the TaxaMatrix along the taxa
        axis. Copy-on-manipulation routine.

        Parameters
        ----------
        values : DensePhasedGenotypeMatrix, numpy.ndarray
            Values to be appended to append to the Matrix.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaMatrix
            A copy of the TaxaMatrix with values appended to the taxa axis
            Note that adjoin does not occur in-place: a new Matrix is allocated
            and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")   # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot adjoin: 'taxa_grp' argument is required")

        # adjoin values
        values = numpy.append(self._mat, values, axis = self.taxa_axis)
        if self._taxa is not None:
            taxa = numpy.append(self._taxa, taxa, axis = 0)
        if self._taxa_grp is not None:
            taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)

        out = self.__class__(
            mat = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            A DenseTaxaMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.delete_taxa(
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
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            A DenseTaxaMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseTaxaMatrix is allocated and filled.
        """
        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # delete values
        mat = numpy.delete(mat, obj, axis = self.taxa_axis)
        if taxa is not None:
            taxa = numpy.delete(taxa, obj, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.delete(taxa_grp, obj, axis = 0)

        out = self.__class__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
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
            If values is a DenseTaxaMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTaxaMatrix
            A DenseTaxaMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaMatrix is allocated and filled.
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
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            A DenseTaxaMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseTaxaMatrix is allocated and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot insert: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")   # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot insert: 'taxa_grp' argument is required")

        # insert values
        values = numpy.insert(self._mat, obj, values, axis = self.taxa_axis)
        if self._taxa is not None:
            taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
        if self._taxa_grp is not None:
            taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)

        # create output
        out = self.__class__(
            mat = values,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    def select(
            self, 
            indices: ArrayLike, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            The output DenseTaxaMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.taxa_axis:
            out = self.select_taxa(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_taxa(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            The output DenseTaxaMatrix with values selected. Note that select does not
            occur in-place: a new DenseTaxaMatrix is allocated and filled.
        """
        # check for array_like
        check_is_array_like(indices, "indices")

        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # select values
        mat = numpy.take(mat, indices, axis = self.taxa_axis)
        if taxa is not None:
            taxa = numpy.take(taxa, indices, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.take(taxa_grp, indices, axis = 0)

        out = self.__class__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    @classmethod
    def concat(
            cls, 
            mats: Sequence, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            The concatenated DenseTaxaMatrix. Note that concat does not occur in-place:
            a new DenseTaxaMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].taxa_axis:
            out = cls.concat_taxa(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_taxa(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseTaxaMatrix':
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
        out : DenseTaxaMatrix
            The concatenated DenseTaxaMatrix. Note that concat does not occur in-place:
            a new DenseTaxaMatrix is allocated and filled.
        """
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # ensure that we have an array_like of length >= 1
        if len(mats) <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DenseTaxaMatrix
        for i,v in enumerate(mats):
            generic_check_isinstance(v, "mats[{0}]".format(i), cls)

        # make sure dimensions are all identical to first element in mats
        if any(m.mat_ndim != mats[0].mat_ndim for m in mats):
            raise ValueError("cannot concat: not all matrices have the same number of dimensions")

        # extract tuple of shapes for testing compatibility
        shape_t = tuple(zip(*[m.mat.shape for m in mats]))

        # test matrix compatibility (same axis length along non-taxa axes)
        for i,v in enumerate(shape_t):                              # for each index,tuple in shape_t
            if (i != mats[0].taxa_axis) and any(l != v[0] for l in v): # if not the taxa axis AND axis lengths are different
                raise ValueError("cannot concat: matrix shapes do not all align along axis {0}".format(i))

        # create matrix lists
        mat_ls = [m.mat for m in mats]
        taxa_ls = [m.taxa for m in mats]
        taxa_grp_ls = [m.taxa_grp for m in mats]

        # process/error check taxa_ls
        if all(e is None for e in taxa_ls):                             # if all elements are None
            taxa_ls = None                                              # replace list with None
        else:                                                           # else at least one matrix does not have a taxa array
            for i,v in enumerate(taxa_ls):                              # for each index,element in taxa_ls
                if v is None:                                           # if element is None
                    ntaxa = shape_t[mats[0].taxa_axis][i]               # get number of taxa
                    taxa_ls[i] = numpy.empty(ntaxa, dtype = "object")   # replace with array of None

        # process/error check taxa_grp_ls
        if all(e is None for e in taxa_grp_ls):         # if all elements are None
            taxa_grp_ls = None                          # replace list with None
        elif any(e is None for e in taxa_grp_ls):       # else if any elements are None
            raise ValueError("cannot concat: 'taxa_grp' needed for all Matrix in list")

        # concatenate mat, taxa, taxa_grp items
        mat = numpy.concatenate(mat_ls, axis = mats[0].taxa_axis)
        taxa = None if taxa_ls is None else numpy.concatenate(taxa_ls, axis = 0)
        taxa_grp = None if taxa_grp_ls is None else numpy.concatenate(taxa_grp_ls, axis = 0)

        # concatenate everything and put into new DenseTaxaMatrix
        out = cls(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

        return out

    ######### Matrix element in-place-manipulation #########
    def append(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseTaxaMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
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
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    def append_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the Matrix along the taxa axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot append: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object") # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot append: 'taxa_grp' argument is required")

        # append values
        self._mat = numpy.append(self._mat, values, axis = self.taxa_axis)
        if self._taxa is not None:
            self._taxa = numpy.append(self._taxa, taxa, axis = 0)
        if self._taxa_grp is not None:
            self._taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)

        # reset metadata
        self._taxa_grp_len = None
        self._taxa_grp_name = None
        self._taxa_grp_stix = None
        self._taxa_grp_spix = None

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
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def remove_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the taxa axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = self.taxa_axis)

        if self._taxa is not None:
            self._taxa = numpy.delete(self._taxa, obj, axis = 0)
        if self._taxa_grp is not None:
            self._taxa_grp = numpy.delete(self._taxa_grp, obj, axis = 0)

        # reset metadata
        self._taxa_grp_len = None
        self._taxa_grp_name = None
        self._taxa_grp_stix = None
        self._taxa_grp_spix = None

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
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
            If values is a DenseTaxaMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
            If values is a DenseTaxaMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis == self.taxa_axis:
            self.incorp(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    def incorp_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the taxa axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        taxa : numpy.ndarray
            Taxa names to incorporate into the Matrix.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if taxa is None:
                taxa = values.taxa
            if taxa_grp is None:
                taxa_grp = values.taxa_grp
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseTaxaMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.taxa_axis) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")  # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot incorp: 'taxa_grp' argument is required")

        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = self.taxa_axis)

        if self._taxa is not None:
            self._taxa = numpy.insert(self._taxa, obj, taxa, axis = 0)
        if self._taxa_grp is not None:
            self._taxa_grp = numpy.insert(self._taxa_grp, obj, taxa_grp, axis = 0)

        # reset metadata
        self._taxa_grp_len = None
        self._taxa_grp_name = None
        self._taxa_grp_stix = None
        self._taxa_grp_spix = None

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
            sort key.
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
        else:
            raise ValueError("cannot lexsort along axis {0}".format(axis))

        return indices

    def lexsort_taxa(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys along the taxa
        axis.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : A (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        # default error message
        emess = "no available keys to sort"

        # if no keys were provided, set a default
        if keys is None:
            keys = (self._taxa, self._taxa_grp)         # taxa default keys
            emess = "taxa, taxa_grp are None"           # taxa error message

        # remove None keys
        keys = tuple(k for k in keys if k is not None)

        # raise error if no keys remain
        if len(keys) == 0:
            raise ValueError("cannot lexsort on axis {0} (taxa axis): {1}".format(self.taxa_axis, emess))

        # raise error if keys are of incompatible length
        if any(len(k) != self.ntaxa for k in keys):
            emess = "keys are not all length {0}".format(self.ntaxa)
            raise ValueError("cannot lexsort on axis {0} (taxa axis): {1}".format(self.taxa_axis, emess))

        # get indices
        indices = numpy.lexsort(keys)

        # return indices
        return indices

    def reorder(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Reorder values in the ``DenseTaxaMatrix``.

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
            self.reorder(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def reorder_taxa(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Reorder elements of the Matrix along the taxa axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : A (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # build a tuple to slice the matrix
        ix = tuple(indices if i == self.taxa_axis else slice(None) for i in range(self.mat_ndim))

        # reorder arrays
        self._mat = self._mat[ix]

        if self._taxa is not None:
            self._taxa = self._taxa[indices]                # reorder taxa array
        if self._taxa_grp is not None:
            self._taxa_grp = self._taxa_grp[indices]        # reorder taxa group array

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
            sort key.
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
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    def sort_taxa(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the Matrix along the taxa axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        # reset taxa group metadata
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

        # get indices for sort
        indices = self.lexsort_taxa(keys, **kwargs)

        # reorder internals
        self.reorder_taxa(indices, **kwargs)

    ################### Grouping Methods ###################
    def group(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Sort the DenseTaxaMatrix along an axis, then populate grouping indices.

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
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

    def group_taxa(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Sort the Matrix along the taxa axis, then populate grouping indices for
        the taxa axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # sort taxa using default keys
        self.sort_taxa()

        if self._taxa_grp is not None:
            # get unique taxa group names, starting indices, group lengths
            uniq = numpy.unique(self._taxa_grp, return_index = True, return_counts = True)
            # make assignments to instance data
            self._taxa_grp_name, self._taxa_grp_stix, self._taxa_grp_len = uniq
            # calculate stop indices
            self._taxa_grp_spix = self._taxa_grp_stix + self._taxa_grp_len

    def ungroup(
            self,
            axis: int = -1,
            **kwargs: dict
        ) -> None:
        """
        Ungroup the DenseTaxaMatrix along an axis by removing grouping metadata.

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
        else:
            raise ValueError("cannot ungroup along axis {0}".format(axis))

    def ungroup_taxa(
            self,
            **kwargs: dict
        ) -> None:
        """
        Ungroup the DenseTaxaMatrix along the taxa axis by removing taxa group 
        metadata.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        # set taxa metadata to None
        self.taxa_grp_name = None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_len = None

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
            The axis along which to test if values are grouped.
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
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped

    def is_grouped_taxa(
            self, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the Matrix has been sorted and grouped along the taxa
        axis.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        return (
            (self._taxa_grp_name is not None) and
            (self._taxa_grp_stix is not None) and
            (self._taxa_grp_spix is not None) and
            (self._taxa_grp_len is not None)
        )

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write ``DenseTaxaMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseTaxaMatrix`` data is stored.
            If ``None``, the ``DenseTaxaMatrix`` is written to the base HDF5 group.

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
            "mat"           : self.mat,
            "taxa"          : self.taxa,
            "taxa_grp"      : self.taxa_grp,
            # metadata
            "taxa_grp_name" : self.taxa_grp_name,
            "taxa_grp_stix" : self.taxa_grp_stix,
            "taxa_grp_spix" : self.taxa_grp_spix,
            "taxa_grp_len"  : self.taxa_grp_len,
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
        ) -> 'DenseTaxaMatrix':
        """
        Read ``DenseTaxaMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseTaxaMatrix`` data is stored.
            If ``None``, ``DenseTaxaMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseTaxaMatrix
            A dense matrix read from file.
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
            # metadata
            "taxa_grp_name"     : None,
            "taxa_grp_stix"     : None,
            "taxa_grp_spix"     : None,
            "taxa_grp_len"      : None,
        }

        ##################################
        ### read mandatory data fields ###

        # read mat array (ndarray dtype = int8)
        data["mat"] = h5py_File_read_ndarray(h5file, groupname + "mat")
        
        #################################
        ### read optional data fields ###

        # read taxa array (ndarray dtype = unicode / object)
        if groupname + "taxa" in h5file:
            data["taxa"] = h5py_File_read_ndarray_utf8(h5file, groupname + "taxa")

        # read taxa_grp array (ndarray dtype = any)
        if groupname + "taxa_grp" in h5file:
            data["taxa_grp"] = h5py_File_read_ndarray(h5file, groupname + "taxa_grp")
        
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
        )

        # copy metadata
        out.taxa_grp_name    = data["taxa_grp_name"]
        out.taxa_grp_stix    = data["taxa_grp_stix"]
        out.taxa_grp_spix    = data["taxa_grp_spix"]
        out.taxa_grp_len     = data["taxa_grp_len"]

        return out



################################## Utilities ###################################
def check_is_DenseTaxaMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseTaxaMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseTaxaMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseTaxaMatrix.__name__,type(v).__name__))
