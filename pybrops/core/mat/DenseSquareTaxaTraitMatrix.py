"""
Module implementing a dense matrix with taxa axes that are square and a trait
axis which is not square, and associated error checking routines.
"""

__all__ = [
    "DenseSquareTaxaTraitMatrix",
    "check_is_DenseSquareTaxaTraitMatrix",
]

import copy
from functools import reduce
from numbers import Integral
from pathlib import Path
from typing import Optional
from typing import Sequence
from typing import Union
import numpy
import h5py
from numpy.typing import ArrayLike
import pandas
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_Sequence_all_type, check_is_Integral, check_is_str
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable
from pybrops.core.error.error_value_pandas import check_pandas_DataFrame_has_column_indices
from pybrops.core.error.error_value_python import check_len
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.SquareTaxaTraitMatrix import SquareTaxaTraitMatrix
from pybrops.core.util.array import flattenix, get_axis
from pybrops.core.util.h5py import h5py_File_write_dict

class DenseSquareTaxaTraitMatrix(
        DenseSquareTaxaMatrix,
        DenseTraitMatrix,
        SquareTaxaTraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
    ):
    """
    A concrete class for dense matrices with taxa axes that are square and a 
    trait axis which is not square.

    The purpose of this abstract class is to merge the following implementations
    and interfaces:

        1. DenseSquareTaxaMatrix (implementation)
        2. DenseTraitMatrix (implementation)
        3. SquareTaxaTraitMatrix (interface)
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
        Constructor for the DenseSquareTaxaTraitMatrix concrete class.

        Parameters
        ----------
        mat : numpy.ndarray
            Matrix used to construct the object.

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

    def __copy__(
            self
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
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
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
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

    ############## Square Metadata Properties ##############
    @DenseSquareTaxaMatrix.square_taxa_axes.getter
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        # square taxa axes are everything except the last axis
        return tuple(range(self._mat.ndim-1))

    #################### Trait metadata ####################
    @DenseTraitMatrix.trait_axis.getter
    def trait_axis(self) -> int:
        """Axis along which traits are stored."""
        # traits are always stored on the last axis.
        return self._mat.ndim - 1

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    def copy(
            self
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Make a shallow copy of the Matrix.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
            A shallow copy of the original DenseSquareTaxaTraitMatrix.
        """
        return copy.copy(self)

    def deepcopy(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Make a deep copy of the Matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
            A deep copy of the original DenseSquareTaxaTraitMatrix.
        """
        return copy.deepcopy(self, memo)

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        taxa : numpy.ndarray
            Taxa names to adjoin to the Matrix.
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray, None
            Trait names to adjoin to the Matrix.
            If ``values`` is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
            A copy of DenseSquareTaxaTraitMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseSquareTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_axes:
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

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
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
        out : DenseSquareTaxaTraitMatrix
            A DenseSquareTaxaTraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseSquareTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_axes:
            out = self.delete_taxa(obj = obj, **kwargs)
        elif axis == self.trait_axis:
            out = self.delete_trait(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None,
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
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
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray, None
            Trait names to insert into the Matrix.
            If ``values`` is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
            A DenseSquareTaxaTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseSquareTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_axes:
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

    def select(
            self, 
            indices: ArrayLike, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
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
        out : DenseSquareTaxaTraitMatrix
            The output DenseSquareTaxaTraitMatrix with values selected. Note that select does not
            occur in-place: a new DenseSquareTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_axes:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis == self.trait_axis:
            out = self.select_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    @classmethod
    def concat(
            cls, 
            mats: Sequence, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
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
        out : DenseSquareTaxaTraitMatrix
            The concatenated DenseSquareTaxaTraitMatrix. Note that concat does not occur in-place:
            a new DenseSquareTaxaTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis in mats[0].square_axes:
            out = cls.concat_taxa(mats = mats, **kwargs)
        elif axis == mats[0].trait_axis:
            out = cls.concat_trait(mats = mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

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
        values : DenseSquareTaxaTraitMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray, None
            Trait names to append to the Matrix.
            If ``values`` is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis in self.square_axes:
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

        if axis in self.square_axes:
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
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
            If values is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray, None
            Trait names to incorporate into the Matrix.
            If ``values`` is a DenseSquareTaxaTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis in self.square_axes:
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
        if axis in self.square_axes:
            indices = self.lexsort_taxa(keys = keys, **kwargs)
        elif axis == self.trait_axis:
            indices = self.lexsort_trait(keys = keys, **kwargs)
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

        if axis in self.square_axes:
            self.reorder_taxa(indices = indices, **kwargs)
        elif axis == self.trait_axis:
            self.reorder_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def sort(
            self, 
            keys: Optional[Union[tuple,numpy.ndarray]] = None, 
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
        if axis in self.square_axes:
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
        Sort the DenseSquareTaxaTraitMatrix along an axis, then populate 
        grouping indices.

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
        if axis in self.square_axes:
            self.group_taxa(**kwargs)
        elif axis == self.trait_axis:
            raise ValueError("cannot group along trait axis {0}".format(axis))
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

    def ungroup(
            self,
            axis: int = -1,
            **kwargs: dict
        ) -> None:
        """
        Ungroup the DenseSquareTaxaMatrix along an axis by removing grouping 
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
        if axis in self.square_axes:
            self.ungroup_taxa(**kwargs)
        elif axis == self.trait_axis:
            raise ValueError("cannot ungroup along trait axis {0}".format(axis))
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
            Axis to test for grouping.
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

        if axis in self.square_axes:
            grouped = self.is_grouped_taxa(**kwargs)
        elif axis == self.trait_axis:
            raise ValueError("cannot test for grouping along trait axis {0}".format(axis))
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
        Write DenseSquareTaxaTraitMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseSquareTaxaTraitMatrix`` data is stored.
            If ``None``, the ``DenseSquareTaxaTraitMatrix`` is written to the base HDF5 group.

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
            "trait"         : self.trait,
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

    # TODO: remove me
    def to_pandas(
            self,
            taxa_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            trait_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            value_colname: str = "value",
            **kwargs: dict
        ) -> pandas.DataFrame:
        """
        Export a DenseSquareTaxaTraitMatrix to a pandas.DataFrame.

        Parameters
        ----------
        taxa_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``None``, then do not export taxa columns.
        
        taxa_grp_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``None``, then do not export taxa group columns.

        trait_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one column name is assumed.
            If ``None``, then do not export trait column.
        
        value_colname : str, default = "value"
            Name of the value column.
        
        kwargs : dict
            Additional keyword arguments to be passed for dataframe construction.

        Returns
        -------
        out : pandas.DataFrame
            Output pandas dataframe.
        """
        ###
        ### Process inputs
        ###

        # get the number of taxa axes
        ntaxaaxes = self.nsquare_taxa

        # get the number of trait axes
        ntraitaxes = 1

        ### process inputs: taxa_colnames

        # process ``None``
        if taxa_colnames is None:
            taxa_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_colnames, bool):
            taxa_colnames = ["taxa_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)] if taxa_colnames else [None for i in range(ntaxaaxes)]
        elif isinstance(taxa_colnames, str):
            taxa_colnames = [taxa_colnames]
        elif isinstance(taxa_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_colnames`` must be of type ``bool``, ``str``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_colnames, "taxa_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_colnames, "taxa_colnames", (str,type(None)))
        
        ### process inputs: taxa_grp_colnames

        # process ``None``
        if taxa_grp_colnames is None:
            taxa_grp_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_grp_colnames, bool):
            taxa_grp_colnames = ["taxa_grp_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)] if taxa_grp_colnames else [None for i in range(ntaxaaxes)]
        elif isinstance(taxa_grp_colnames, str):
            taxa_grp_colnames = [taxa_grp_colnames]
        elif isinstance(taxa_grp_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_grp_colnames`` must be of type ``bool``, ``str``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_grp_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_grp_colnames, "taxa_grp_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_grp_colnames, "taxa_grp_colnames", (str,type(None)))

        ### process inputs: trait_colnames

        # process ``None``
        if trait_colnames is None:
            trait_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(trait_colnames, bool):
            trait_colnames = ["trait_"+str(i).zfill(len(str(ntraitaxes))) for i in range(ntraitaxes)] if trait_colnames else [None for i in range(ntraitaxes)]
        elif isinstance(trait_colnames, str):
            trait_colnames = [trait_colnames]
        elif isinstance(trait_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``trait_colnames`` must be of type ``bool``, ``str``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(trait_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(trait_colnames, "trait_colnames", ntraitaxes)

        # check sequence element types
        check_Sequence_all_type(trait_colnames, "trait_colnames", (str,type(None)))

        ### process inputs: value_colname
        check_is_str(value_colname, "value_colname")

        ###
        ### Process data
        ###

        # get taxa names
        taxa = numpy.array(["Taxon"+str(e).zfill(len(str(self.ntaxa))) for e in range(self.ntaxa)], dtype = object) if self.taxa is None else self.taxa

        # get taxa_grp names
        taxa_grp = numpy.array([pandas.NA for i in range(self.ntaxa)], dtype=object) if self.taxa_grp is None else self.taxa_grp

        # get trait names
        trait = numpy.array(["Trait"+str(i).zfill(len(str(self.ntrait))) for i in range(self.ntrait)], dtype=object) if self.trait is None else self.trait

        print(self.ntrait)

        # calculate flattened array and corresponding axis indices
        flatmat, axisix_data = flattenix(self.mat)

        ###
        ### Export data to dict then to pandas and return
        ###

        # create empty dict
        out_dict = {}

        # save taxa column data
        for ix,taxa_colname in zip(self.square_taxa_axes,taxa_colnames):
            if taxa_colname is None:
                continue
            axisix = axisix_data[ix]
            out_dict.update({taxa_colname: taxa[axisix]})

        # save taxa grp column data
        for ix,taxa_grp_colname in zip(self.square_taxa_axes,taxa_grp_colnames):
            if taxa_grp_colname is None:
                continue
            axisix = axisix_data[ix]
            out_dict.update({taxa_grp_colname: taxa_grp[axisix]})
        
        # save trait data
        trait_colname = trait_colnames[0]
        if trait_colname is not None:
            axisix = axisix_data[self.trait_axis]
            out_dict.update({trait_colname: trait[axisix]})
        
        # save values
        out_dict.update({value_colname: flatmat})

        # create a pandas DataFrame from the data
        out = pandas.DataFrame(out_dict, **kwargs)

        return out

    # TODO: remove me
    def to_csv(
            self,
            filename: str,
            taxa_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            trait_colnames: Optional[Union[bool,str,Sequence[Union[str,None]]]] = True,
            value_colname: str = "value",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Export a DenseSquareTaxaTraitMatrix to a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name to which to write.
        
        taxa_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``None``, then do not export taxa columns.
        
        taxa_grp_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``None``, then do not export taxa group columns.

        trait_colnames : bool, str, Sequence of str or None, None, default = True
            Sequence of column names for which to name trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then assign default names if ``True`` or do not export columns if ``False``.
            If ``str``, then one column name is assumed.
            If ``None``, then do not export trait column.
        
        value_colname : str, default = "value"
            Name of the value column.
        
        sep : str, default = ","
            Separator to use in the exported CSV file.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV file.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # convert DenseTwoWayDHAdditiveProgenyGeneticCovarianceMatrix to pandas.DataFrame
        df = self.to_pandas(
            taxa_colnames = taxa_colnames,
            taxa_grp_colnames = taxa_grp_colnames,
            trait_colnames = trait_colnames,
            value_colname = value_colname,
        )

        # export using pandas
        df.to_csv(
            path_or_buf = filename,
            sep         = sep,
            header      = header,
            index       = index,
            **kwargs
        )

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Read DenseSquareTaxaTraitMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which DenseSquareTaxaTraitMatrix data is stored.
            If None, DenseSquareTaxaTraitMatrix is read from base HDF5 group.

        Returns
        -------
        out : DenseSquareTaxaTraitMatrix
            A dense matrix read from file.
        """
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
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an h5py.File, make sure mode is in at least ``r`` mode, and copy pointer
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
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = ["mat"]                               # all required arguments

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)

        ########################################################
        ### read data from HDF5 file and (optionally) close ####
        
        # output dictionary
        data = {
            "mat"           : None,
            "taxa"          : None,
            "taxa_grp"      : None,
            "trait"         : None,
            # metadata
            "taxa_grp_name" : None,
            "taxa_grp_stix" : None,
            "taxa_grp_spix" : None,
            "taxa_grp_len"  : None,
        }
        
        ##################################
        ### read mandatory data fields ###

        # read mat array
        data["mat"] = h5file[groupname + "mat"][()]
        
        #################################
        ### read optional data fields ###

        # read taxa data, if "groupname/taxa" in HDF5 file
        if groupname + "taxa" in h5file:
            data["taxa"] = numpy.array([s.decode("utf-8") if isinstance(s,bytes) else s for s in h5file[groupname+"taxa"][()]], dtype=object)

        # read taxa_grp data, if "groupname/taxa_grp" in HDF5 file
        if groupname + "taxa_grp" in h5file:
            data["taxa_grp"] = h5file[groupname + "taxa_grp"][()]
        
        # read trait data, if "groupname/trait" in HDF5 file
        if groupname + "trait" in h5file:
            data["trait"] = numpy.array([s.decode("utf-8") if isinstance(s,bytes) else s for s in h5file[groupname+"trait"][()]], dtype=object)

        #####################################
        ### read optional metadata fields ###

        # read taxa_grp_name data, if "groupname/taxa_grp_name" in HDF5 file
        if groupname + "taxa_grp_name" in h5file:
            data["taxa_grp_name"] = h5file[groupname + "taxa_grp_name"][()]

        # read taxa_grp_stix data, if "groupname/taxa_grp_stix" in HDF5 file
        if groupname + "taxa_grp_stix" in h5file:
            data["taxa_grp_stix"] = h5file[groupname + "taxa_grp_stix"][()]

        # read taxa_grp_spix data, if "groupname/taxa_grp_spix" in HDF5 file
        if groupname + "taxa_grp_spix" in h5file:
            data["taxa_grp_spix"] = h5file[groupname + "taxa_grp_spix"][()]

        # read taxa_grp_len data, if "groupname/taxa_grp_len" in HDF5 file
        if groupname + "taxa_grp_len" in h5file:
            data["taxa_grp_len"] = h5file[groupname + "taxa_grp_len"][()]

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ################### Object creation ####################
        
        # create object from read data
        mat = cls(
            mat      = data["mat"],
            taxa     = data["taxa"],
            taxa_grp = data["taxa_grp"],
            trait    = data["trait"],
        )

        # assign metadata
        mat.taxa_grp_name = data["taxa_grp_name"]
        mat.taxa_grp_stix = data["taxa_grp_stix"]
        mat.taxa_grp_spix = data["taxa_grp_spix"]
        mat.taxa_grp_len  = data["taxa_grp_len"]

        return mat

    # TODO: remove me
    @classmethod
    def from_pandas(
            cls,
            df: pandas.DataFrame,
            taxa_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            trait_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            value_colname: Union[str,Integral] = "value",
            ntaxaaxes: Integral = 2,
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Import a DenseSquareTaxaTraitMatrix from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.

        taxa_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read taxa axis columns.
            There must be one column name or index for each taxa axis.
            If an element in the sequence cannot be ``None``; all columns must be present.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``Integral``, then one taxa axis column is assumed.
            If ``None``, then raise error.
        
        taxa_grp_colnames : bool, str, Integral, Sequence of str or None, None, default = True
            Sequence of column names from which to read taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``Integral``, then one taxa group axis column is assumed.
            If ``None``, then raise error.

        trait_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one column name is assumed.
            If ``Integral``, then one column name is assumed.
            If ``None``, then raise error.
        
        value_colname : str, Integral, default = "value"
            Name or index of the value column.
        
        ntaxaaxes : Integral
            Expected number of taxa axes for the matrix.
        
        kwargs : dict
            Additional keyword arguments to be passed for dataframe construction.

        Returns
        -------
        out : pandas.DataFrame
            Output pandas dataframe.
        """
        ###
        ### Type checks and preprocessing
        ###

        # get the number of taxa axes
        check_is_Integral(ntaxaaxes, "ntaxaaxes")

        # get the number of trait axes
        ntraitaxes = 1

        ### process inputs: df
        check_is_pandas_DataFrame(df, "df")

        ### process inputs: taxa_colnames

        # process ``None``
        if taxa_colnames is None:
            taxa_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_colnames, bool):
            if taxa_colnames:
                taxa_colnames = ["taxa_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)]
            else:
                raise ValueError("must provide taxa column name or index values")
        elif isinstance(taxa_colnames, (str,Integral)):
            taxa_colnames = [taxa_colnames]
        elif isinstance(taxa_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_colnames`` must be of type ``bool``, ``str``, ``Integral``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_colnames, "taxa_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_colnames, "taxa_colnames", (str,Integral))

        # convert taxa_colnames to list so it is mutable
        taxa_colnames = list(taxa_colnames)

        # convert taxa_colnames sequence to all Integral
        for i,elem in enumerate(taxa_colnames):
            if isinstance(elem, str):
                taxa_colnames[i] = df.columns.get_loc(elem)

        # check column indices
        check_pandas_DataFrame_has_column_indices(df, "df", *taxa_colnames)

        ### process inputs: taxa_grp_colnames

        # process ``None``
        if taxa_grp_colnames is None:
            taxa_grp_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(taxa_grp_colnames, bool):
            taxa_grp_colnames = ["taxa_grp_"+str(i).zfill(len(str(ntaxaaxes))) for i in range(ntaxaaxes)] if taxa_grp_colnames else [None for i in range(ntaxaaxes)]
        elif isinstance(taxa_grp_colnames, (str,Integral)):
            taxa_grp_colnames = [taxa_grp_colnames]
        elif isinstance(taxa_grp_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``taxa_grp_colnames`` must be of type ``bool``, ``str``, ``Integral``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(taxa_grp_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(taxa_grp_colnames, "taxa_grp_colnames", ntaxaaxes)

        # check sequence element types
        check_Sequence_all_type(taxa_grp_colnames, "taxa_grp_colnames", (str,Integral,type(None)))

        # convert taxa_grp_colnames to list so it is mutable
        taxa_grp_colnames = list(taxa_grp_colnames)

        # convert taxa_grp_colnames sequence to all Integral
        for i,elem in enumerate(taxa_grp_colnames):
            if isinstance(elem, str):
                taxa_grp_colnames[i] = df.columns.get_loc(elem)

        # create a reduced set that does not contain None (skip None columns)
        reduced = [e for e in taxa_grp_colnames if e is not None]

        # check column indices
        check_pandas_DataFrame_has_column_indices(df, "df", *reduced)

        ### process inputs: trait_colnames

        # process ``None``
        if trait_colnames is None:
            trait_colnames = False
        
        # process ``bool``, ``str``, ``Sequence``
        if isinstance(trait_colnames, bool):
            if trait_colnames:
                trait_colnames = ["trait_"+str(i).zfill(len(str(ntraitaxes))) for i in range(ntraitaxes)]
            else:
                raise ValueError("must provide trait column name or index values")
        elif isinstance(trait_colnames, (str,Integral)):
            trait_colnames = [trait_colnames]
        elif isinstance(trait_colnames, Sequence):
            pass
        else:
            raise TypeError(
                "``trait_colnames`` must be of type ``bool``, ``str``, ``Integral``, ``Sequence``, or ``None`` but received type ``{0}``".format(
                    type(trait_colnames).__name__
                )
            )
        
        # check the sequence length
        check_len(trait_colnames, "trait_colnames", ntraitaxes)

        # check sequence element types
        check_Sequence_all_type(trait_colnames, "trait_colnames", (str,type(None)))

        # convert trait_colnames to list so it is mutable
        trait_colnames = list(trait_colnames)

        # convert trait_colnames sequence to all Integral
        for i,elem in enumerate(trait_colnames):
            if isinstance(elem, str):
                trait_colnames[i] = df.columns.get_loc(elem)

        # check column indices
        check_pandas_DataFrame_has_column_indices(df, "df", *trait_colnames)

        ### process inputs: value_colname
        if isinstance(value_colname, str):
            value_colname = df.columns.get_loc(value_colname)
        
        check_is_Integral(value_colname, "value_colname")

        #####
        ##### Processing
        #####

        ###
        ### Process taxa
        ###

        # get taxa columns data
        taxa_data_ls = [df.iloc[:,ix].to_numpy(dtype=object) for ix in taxa_colnames]

        # get unique taxa and corresponding indices
        taxa_vector_ls, taxa_index_ls = tuple(zip(*[numpy.unique(e, return_index=True) for e in taxa_data_ls]))

        # get unique taxa for all columns: union all taxa together
        taxa = reduce(numpy.union1d, taxa_vector_ls)

        # allocate index arrays for each taxa column
        taxaix_ls = [numpy.empty(len(e), dtype=int) for e in taxa_data_ls]

        # calculate taxa indices
        for i,taxon in enumerate(taxa):
            for taxaix,taxa_data in zip(taxaix_ls,taxa_data_ls):
                taxaix[taxa_data == taxon] = i
        
        ###
        ### Process taxa_grp
        ###

        # get taxa_grp columns data
        taxa_grp_data_ls = [None if ix is None else df.iloc[:,ix].to_numpy(dtype=int) for ix in taxa_grp_colnames]

        # get optional taxa group data
        taxa_grp = None
        for taxa_vector,taxa_index,taxa_grp_data in zip(taxa_vector_ls,taxa_index_ls,taxa_grp_data_ls):
            if taxa_grp_data is not None:
                if taxa_grp is None:
                    taxa_grp = numpy.empty(len(taxa), dtype=int)
                for i, taxon in enumerate(taxa):
                    if taxon in taxa_vector:
                        taxa_grp[i] = taxa_grp_data[(taxa_index[taxa_vector == taxon][0])]

        ###
        ### Process trait
        ###

        # get trait column data
        trait_data = df.iloc[:,(trait_colnames[0])].to_numpy(dtype=object)

        # calculate unique 
        trait, traitix = numpy.unique(trait_data, return_inverse=True)

        # combine trait names
        # trait = reduce(numpy.union1d, blah)

        ###
        ### Process values
        ###

        # get value column data
        value_data = df.iloc[:,(value_colname)].to_numpy(dtype=float)

        # get array dimensions
        ntaxa_tup = tuple(len(taxa) for _ in range(ntaxaaxes))
        ntrait_tup = (len(trait),)
        dim_tup = ntaxa_tup + ntrait_tup

        # allocate NaN array for matrix
        mat = numpy.full(dim_tup, numpy.nan, dtype = float)

        # construct index tuple
        ix_tup = tuple(taxaix_ls) + (traitix,)

        # overwirte NaN values with values
        mat[ix_tup] = value_data

        # construct an object
        out = cls(
            mat = mat, 
            taxa = taxa, 
            taxa_grp = taxa_grp, 
            trait = trait, 
        )

        return out

    # TODO: remove me
    @classmethod
    def from_csv(
            cls,
            filename: str, 
            taxa_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            taxa_grp_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            trait_colnames: Optional[Union[bool,str,Integral,Sequence[Union[str,Integral]]]] = True,
            value_colname: Union[str,Integral] = "value",
            ntaxaaxes: Integral = 2,
            sep: str = ',',
            header: int = 0,
            **kwargs: dict
        ) -> 'DenseSquareTaxaTraitMatrix':
        """
        Read a DenseSquareTaxaTraitMatrix from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.

        taxa_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read taxa axis columns.
            There must be one column name or index for each taxa axis.
            If an element in the sequence cannot be ``None``; all columns must be present.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa axis column is assumed.
            If ``Integral``, then one taxa axis column is assumed.
            If ``None``, then raise error.
        
        taxa_grp_colnames : bool, str, Integral, Sequence of str or None, None, default = True
            Sequence of column names from which to read taxa group axis columns.
            There must be one column name for each taxa axis.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one taxa group axis column is assumed.
            If ``Integral``, then one taxa group axis column is assumed.
            If ``None``, then raise error.

        trait_colnames : bool, str, Integral, Sequence of str or Integral, None, default = True
            Sequence of column names or indices from which to read trait columns.
            There must be one column name for each trait.
            If an element in the sequence is ``None``, then do not export the column.
            If ``bool``, then import using default names if ``True`` or raise error if ``False``.
            If ``str``, then one column name is assumed.
            If ``Integral``, then one column name is assumed.
            If ``None``, then raise error.
        
        value_colname : str, Integral, default = "value"
            Name or index of the value column.
        
        ntaxaaxes : Integral
            Expected number of taxa axes for the matrix.

        sep : str
            CSV file separator.
        
        header : int
            Header row index.
        
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : CSVInputOutput
            An object read from a CSV file.
        """
        # read dataframe from file to pandas
        df = pandas.read_csv(
            filepath_or_buffer = filename,
            sep = sep,
            header = header,
            **kwargs
        )

        # convert pandas to matrix
        out = cls.from_pandas(
            df = df,
            taxa_colnames = taxa_colnames,
            taxa_grp_colnames = taxa_grp_colnames,
            trait_colnames = trait_colnames,
            value_colname = value_colname,
            ntaxaaxes = ntaxaaxes,
        )

        return out

################################## Utilities ###################################
def check_is_DenseSquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSquareTaxaTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSquareTaxaTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseSquareTaxaTraitMatrix.__name__,type(v).__name__))
