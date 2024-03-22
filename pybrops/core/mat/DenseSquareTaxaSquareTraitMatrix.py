"""
Module implementing a dense matrix with taxa and trait axes are square and 
associated error checking routines.
"""

__all__ = [
    "DenseSquareTaxaSquareTraitMatrix",
    "check_is_DenseSquareTaxaSquareTraitMatrix",
]

import copy
from pathlib import Path
from typing import Optional
from typing import Sequence
from typing import Union
import numpy
import h5py
from numpy.typing import ArrayLike
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group
from pybrops.core.error.error_value_h5py import check_h5py_File_is_readable
from pybrops.core.error.error_value_h5py import check_h5py_File_is_writable
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.core.mat.DenseSquareTraitMatrix import DenseSquareTraitMatrix
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.SquareTaxaSquareTraitMatrix import SquareTaxaSquareTraitMatrix
from pybrops.core.mat.util import get_axis
from pybrops.core.util.h5py import h5py_File_read_ndarray
from pybrops.core.util.h5py import h5py_File_read_ndarray_utf8
from pybrops.core.util.h5py import h5py_File_write_dict

class DenseSquareTaxaSquareTraitMatrix(
        DenseSquareTaxaMatrix,
        DenseSquareTraitMatrix,
        SquareTaxaSquareTraitMatrix,
    ):
    """
    docstring for DenseSquareTaxaSquareTraitMatrix.
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
        Constructor for DenseSquareTaxaSquareTraitMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            A numpy.ndarray used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        # store data
        self.mat = mat
        self.taxa = taxa
        self.taxa_grp = taxa_grp
        self.trait = trait

        # set metadate to None
        self.taxa_grp_stix = None
        self.taxa_grp_spix = None
        self.taxa_grp_name = None
        self.taxa_grp_len  = None

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseSquareTaxaSquareTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseSquareTaxaSquareTraitMatrix
            A shallow copy of the matrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            trait = copy.copy(self.trait),
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
        ) -> 'DenseSquareTaxaSquareTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquareTaxaSquareTraitMatrix
            A deep copy of the matrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            trait = copy.deepcopy(self.trait, memo),
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        return out

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    # mat                       (inherited from DenseSquareMatrix)
    # mat_ndim                  (inherited from DenseSquareMatrix)
    # mat_shape                 (inherited from DenseSquareMatrix)

    ############## Square Metadata Properties ##############
    @property
    def nsquare(self) -> int:
        """Number of axes that are square."""
        return self.nsquare_taxa + self.nsquare_trait

    @property
    def square_axes(self) -> tuple:
        """Axis indices for axes that are square."""
        return self.square_taxa_axes + self.square_trait_axes

    @property
    def square_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        return self.square_taxa_axes_len + self.square_trait_axes_len

    # nsquare_taxa              (inherited from DenseSquareTaxaMatrix)
    
    @property
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        return (0,1)
    
    # square_taxa_axes_len      (inherited from DenseSquareTaxaMatrix)

    # nsquare_trait             (inherited from DenseSquareTraitMatrix)
    
    @property
    def square_trait_axes(self) -> tuple:
        """Axis indices for trait axes that are square."""
        return (2,3)
    
    # square_trait_axes_len     (inherited from DenseSquareTraitMatrix)

    ################# Taxa Data Properites #################

    ############### Taxa Metadata Properites ###############

    ################### Fill data lookup ###################
    # _fill_value               (inherited from DenseSquareMatrix)

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    def is_square(
            self
        ) -> bool:
        """
        Determine whether the axis lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square axes are the same length.
            ``False`` if not all square axes are the same length.
        """
        return self.is_square_taxa() and self.is_square_trait()

    # is_square_taxa            (inherited from DenseSquareTaxaMatrix)
    # is_square_trait           (inherited from DenseSquareTraitMatrix)

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
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
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A copy of DenseSquareTraitMatrix with values adjoined to axis. Note that adjoin does
            not occur in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_taxa_axes:
            out = self.adjoin_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis in self.square_trait_axes:
            out = self.adjoin_trait(
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot adjoin along axis {0}".format(axis))

        return out

    # adjoin_taxa               (inherited from DenseSquareTaxaMatrix)
    # adjoin_trait              (inherited from DenseSquareTraitMatrix)

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
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
        out : DenseSquareTraitMatrix
            A DenseSquareTraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_taxa_axes:
            out = self.delete_taxa(obj = obj, **kwargs)
        elif axis in self.square_trait_axes:
            out = self.delete_trait(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    # delete_taxa               (inherited from DenseSquareTaxaMatrix)
    # delete_trait              (inherited from DenseSquareTraitMatrix)

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
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
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to insert into the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTaxaSquareTraitMatrix
            A DenseSquareTaxaSquareTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseSquareTaxaSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_taxa_axes:
            out = self.insert_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis in self.square_trait_axes:
            out = self.insert_trait(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    # insert_taxa               (inherited from DenseSquareTaxaMatrix)
    # insert_trait              (inherited from DenseSquareTraitMatrix)

    def select(
            self, 
            indices: ArrayLike, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
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
        out : DenseSquareTraitMatrix
            The output DenseSquareTraitMatrix with values selected. Note that select does not
            occur in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_taxa_axes:
            out = self.select_taxa(indices = indices, **kwargs)
        elif axis in self.square_trait_axes:
            out = self.select_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    # select_taxa               (inherited from DenseSquareTaxaMatrix)
    # select_trait              (inherited from DenseSquareTraitMatrix)

    @classmethod
    def concat(
            cls, 
            mats: Sequence, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
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
        out : DenseSquareTraitMatrix
            The concatenated DenseSquareTraitMatrix. Note that concat does not occur in-place:
            a new DenseSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis in mats[0].square_taxa_axes:
            out = cls.concat_taxa(mats, **kwargs)
        elif axis in mats[0].square_trait_axes:
            out = cls.concat_trait(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    # concat_taxa               (inherited from DenseSquareTaxaMatrix)
    # concat_trait              (inherited from DenseSquareTraitMatrix)

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
        values : DenseHaplotypeMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        taxa : numpy.ndarray
            Taxa names to append to the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to append to the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to append to the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis in self.square_taxa_axes:
            self.append_taxa(
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis in self.square_trait_axes:
            self.append_trait(
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    # append_taxa               (inherited from DenseSquareTaxaMatrix)
    # append_trait              (inherited from DenseSquareTraitMatrix)

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

        # dispatch functions
        if axis in self.square_taxa_axes:
            self.remove_taxa(obj = obj, **kwargs)
        elif axis in self.square_trait_axes:
            self.remove_trait(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    # remove_taxa               (inherited from DenseSquareTaxaMatrix)
    # remove_trait              (inherited from DenseSquareTraitMatrix)

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
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to incorporate into the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        trait : numpy.ndarray
            Trait names to incorporate into the Matrix.
            If values is a DenseSquareTaxaSquareTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch function
        if axis in self.square_taxa_axes:
            self.incorp_taxa(
                obj = obj,
                values = values,
                taxa = taxa,
                taxa_grp = taxa_grp,
                **kwargs
            )
        elif axis in self.square_trait_axes:
            self.incorp_trait(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    # incorp_taxa               (inherited from DenseSquareTaxaMatrix)
    # incorp_trait              (inherited from DenseSquareTraitMatrix)

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
        if axis in self.square_taxa_axes:
            indices = self.lexsort_taxa(keys = keys, **kwargs)
        elif axis in self.square_trait_axes:
            indices = self.lexsort_trait(keys = keys, **kwargs)
        else:
            raise ValueError("cannot lexsort along axis {0}".format(axis))

        return indices

    # lexsort_taxa              (inherited from DenseSquareTaxaMatrix)
    # lexsort_trait             (inherited from DenseSquareTraitMatrix)

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

        if axis in self.square_taxa_axes:
            self.reorder_taxa(indices = indices, **kwargs)
        elif axis in self.square_trait_axes:
            self.reorder_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    # reorder_taxa              (inherited from DenseSquareTaxaMatrix)
    # reorder_trait             (inherited from DenseSquareTraitMatrix)

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
        if axis in self.square_taxa_axes:
            self.sort_taxa(keys = keys, **kwargs)
        elif axis in self.square_trait_axes:
            self.sort_trait(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    # sort_taxa                 (inherited from DenseSquareTaxaMatrix)
    # sort_trait                (inherited from DenseSquareTraitMatrix)

    ################### Grouping Methods ###################
    def group(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Sort matrix along axis, then populate grouping indices for the axis.

        Parameters
        ----------
        axis : int
            Axis along which values should be grouped.
        kwargs : dict
            Additional keyword arguments.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis in self.square_taxa_axes:
            self.group_taxa(**kwargs)
        elif axis in self.square_trait_axes:
            raise ValueError("cannot group along trait axis {0}: trait axes are not groupable".format(axis))
        else:
            raise ValueError("cannot group along axis {0}".format(axis))

    # group_taxa is unaltered

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
        if axis in self.square_taxa_axes:
            self.ungroup_taxa(**kwargs)
        elif axis in self.square_trait_axes:
            raise ValueError("cannot ungroup along trait axis {0}: trait axes are not groupable".format(axis))
        else:
            raise ValueError("cannot ungroup along axis {0}".format(axis))

    # ungroup_taxa is unaltered

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
            Axis along which to test if values are grouped.
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

        if axis in self.square_taxa_axes:
            grouped = self.is_grouped_taxa(**kwargs)
        elif axis in self.square_trait_axes:
            raise ValueError("cannot test for grouping along trait axis {0}: trait axes are not groupable".format(axis))
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped

    # is_grouped_taxa is unaltered

    ################### Matrix File I/O ####################
    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write ``DenseSquareTaxaSquareTraitMatrix`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which the ``DenseSquareTaxaSquareTraitMatrix`` data is stored.
            If ``None``, the ``DenseSquareTaxaSquareTraitMatrix`` is written to the base HDF5 group.

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

    ############################## Class Methods ###############################

    ################### Matrix File I/O ####################
    @classmethod
    def from_hdf5(
            cls, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseSquareTaxaSquareTraitMatrix':
        """
        Read ``DenseSquareTaxaSquareTraitMatrix`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseSquareTaxaSquareTraitMatrix`` data is stored.
            If ``None``, ``DenseSquareTaxaSquareTraitMatrix`` is read from base HDF5 group.

        Returns
        -------
        out : DenseSquareTaxaSquareTraitMatrix
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
        
        # read trait array (ndarray dtype = unicode / object)
        if groupname + "trait" in h5file:
            data["trait"] = h5py_File_read_ndarray_utf8(h5file, groupname + "trait")

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
            trait       = data["trait"],
        )

        # copy metadata
        out.taxa_grp_name   = data["taxa_grp_name"]
        out.taxa_grp_stix   = data["taxa_grp_stix"]
        out.taxa_grp_spix   = data["taxa_grp_spix"]
        out.taxa_grp_len    = data["taxa_grp_len"]

        return out



################################## Utilities ###################################
def check_is_DenseSquareTaxaSquareTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSquareTaxaSquareTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSquareTaxaSquareTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseSquareTaxaSquareTraitMatrix.__name__,type(v).__name__))
