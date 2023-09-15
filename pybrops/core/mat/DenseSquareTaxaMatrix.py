"""
Module implementing a dense matrix with taxa metadata and axes that are square
and associated error checking routines.
"""

__all__ = [
    "DenseSquareTaxaMatrix",
    "check_is_DenseSquareTaxaMatrix",
]

from typing import Optional, Sequence, Union
import numpy
from numpy.typing import ArrayLike

from pybrops.core.error.error_type_python import check_is_array_like
from pybrops.core.error.error_attr_python import check_is_iterable
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.error.error_generic_python import generic_check_isinstance
from pybrops.core.mat.DenseSquareMatrix import DenseSquareMatrix
from pybrops.core.mat.DenseTaxaMatrix import DenseTaxaMatrix
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix

class DenseSquareTaxaMatrix(DenseSquareMatrix,DenseTaxaMatrix,SquareTaxaMatrix):
    """
    A concrete class for dense matrices with taxa metadata and axes that are
    square.

    The purpose of this abstract class is to merge the following implementations
    and interfaces:

        1. DenseSquareMatrix (implementation)
        2. DenseTaxaMatrix (implementation)
        3. SquareTaxaMatrix (interface)
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
        Constructor for the DenseSquareTaxaMatrix concrete class.

        Parameters
        ----------
        mat : numpy.ndarray
            Matrix used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseSquareTaxaMatrix, self).__init__(
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

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    # mat inherited from DenseSquareMatrix

    ############## Square Metadata Properties ##############
    # nsquare inherited from DenseSquareMatrix

    # square_axes inherited from DenseSquareMatrix

    # square_axes_len inherited from DenseSquareMatrix

    ################# Taxa Data Properites #################

    ############### Taxa Metadata Properites ###############

    ################### Fill data lookup ###################
    # _fill_value inherited from DenseSquareMatrix

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    # is_square inherited from DenseSquareMatrix

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaMatrix':
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTaxaMatrix
            A copy of DenseSquareTaxaMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseSquareTaxaMatrix is allocated and filled.
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
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_taxa(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

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
        out : DenseSquareTaxaMatrix
            A copy of DenseSquareTaxaMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseSquareTaxaMatrix is allocated and filled.
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
            if (i not in self.square_axes) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")   # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot adjoin: 'taxa_grp' argument is required")

        # calculate shape of new matrix to allocate
        mshape = list(self.mat_shape)   # get matrix shape of self
        vshape = values.shape           # get matrix shape of values
        for ix in self.square_axes:     # for each square axis index
            mshape[ix] += vshape[ix]    # add value axis length to self shape
        mshape = tuple(mshape)          # convert list to tuple

        # allocate memory and fill with default fill values
        mat = numpy.empty(mshape, dtype = self._mat.dtype)
        mat.fill(self._fill_value[mat.dtype])

        # calculate indexing tuples for copying matrix values into 'mat'
        nd = mat.ndim           # number of dimensions
        sms = self.mat_shape    # matrix shape for self
        ssa = self.square_axes  # tuple of square taxa axes
        ms = mat.shape          # output matrix shape

        # calculate indexing tuple for self._mat
        # if axis is a square axis then
        #     slice to select 0 to 'a' along axis (numpy: [...,0:a,...])
        # else
        #     slice everything (numpy: [...,:,...])
        six = tuple(slice(0,sms[i]) if i in ssa else slice(None) for i in range(nd))

        # calculate indexing tuple for self._mat
        # if axis is a square axis then
        #     slice to select 'a' to 'b' along axis (numpy: [...,a:b,...])
        # else
        #     slice everything (numpy: [...,:,...])
        vix = tuple(slice(sms[i],ms[i]) if i in ssa else slice(None) for i in range(nd))

        # adjoin values by copying
        mat[six] = self._mat    # copy self into matrix
        mat[vix] = values       # copy values into matrix

        # adjoin values
        if self._taxa is not None:
            taxa = numpy.append(self._taxa, taxa, axis = 0)
        if self._taxa_grp is not None:
            taxa_grp = numpy.append(self._taxa_grp, taxa_grp, axis = 0)

        # create output
        out = self.__class__(
            mat = mat,
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
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            A DenseSquareTaxaMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseSquareTaxaMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_axes:
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
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            A DenseSquareTaxaMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseSquareTaxaMatrix is allocated and filled.
        """
        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # delete values
        for axis in self.square_axes:
            mat = numpy.delete(mat, obj, axis = axis)
        if taxa is not None:
            taxa = numpy.delete(taxa, obj, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.delete(taxa_grp, obj, axis = 0)

        # create output
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
        ) -> 'DenseSquareTaxaMatrix':
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
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa field, providing this argument overwrites the field.
        taxa_grp : numpy.ndarray
            Taxa groups to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            taxa_grp field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTaxaMatrix
            A DenseSquareTaxaMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseSquareTaxaMatrix is allocated and filled.
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
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    # TODO: FIXME
    def insert_taxa(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            A DenseSquareTaxaMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseSquareTaxaMatrix is allocated and filled.
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
            if (i not in self.square_axes) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")   # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot insert: 'taxa_grp' argument is required")

        # TODO: # FIXME: figure out insertion logic
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
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            The output DenseSquareTaxaMatrix with values selected. Note that select does not
            occur in-place: a new DenseSquareTaxaMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_axes:
            out = self.select_taxa(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_taxa(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            The output DenseSquareTaxaMatrix with values selected. Note that select does not
            occur in-place: a new DenseSquareTaxaMatrix is allocated and filled.
        """
        # check for array_like
        check_is_array_like(indices, "indices")

        # get values
        mat = self._mat
        taxa = self._taxa
        taxa_grp = self._taxa_grp

        # select values
        for axis in self.square_axes:
            mat = numpy.take(mat, indices, axis = axis)
        if taxa is not None:
            taxa = numpy.take(taxa, indices, axis = 0)
        if taxa_grp is not None:
            taxa_grp = numpy.take(taxa_grp, indices, axis = 0)

        # create output
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
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            The concatenated DenseSquareTaxaMatrix. Note that concat does not occur in-place:
            a new DenseSquareTaxaMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis in mats[0].square_axes:
            out = cls.concat_taxa(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    # TODO: FIXME
    @classmethod
    def concat_taxa(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseSquareTaxaMatrix':
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
        out : DenseSquareTaxaMatrix
            The concatenated DenseSquareTaxaMatrix. Note that concat does not occur in-place:
            a new DenseSquareTaxaMatrix is allocated and filled.
        """
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # ensure that we have an array_like of length >= 1
        if len(mats) <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DenseHaplotypeMatrix
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

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseHaplotypeMatrix
        # use first element as source of variant data
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
        values : DenseHaplotypeMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
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
            if (i not in self.square_axes) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object") # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot append: 'taxa_grp' argument is required")

        # calculate shape of new matrix to allocate
        mshape = list(self.mat_shape)   # get matrix shape of self
        vshape = values.shape           # get matrix shape of values
        for ix in self.square_axes:     # for each square axis index
            mshape[ix] += vshape[ix]    # add value axis length to self shape
        mshape = tuple(mshape)          # convert list to tuple

        # allocate memory and fill with default fill values
        mat = numpy.empty(mshape, dtype = self._mat.dtype)
        mat.fill(self._fill_value[mat.dtype])

        # calculate indexing tuples for copying matrix values into 'mat'
        nd = mat.ndim           # number of dimensions
        sms = self.mat_shape    # matrix shape for self
        ssa = self.square_axes  # tuple of square taxa axes
        ms = mat.shape          # output matrix shape

        # calculate indexing tuple for self._mat
        # if axis is a square axis then
        #     slice to select 0 to 'a' along axis (numpy: [...,0:a,...])
        # else
        #     slice everything (numpy: [...,:,...])
        six = tuple(slice(0,sms[i]) if i in ssa else slice(None) for i in range(nd))

        # calculate indexing tuple for self._mat
        # if axis is a square axis then
        #     slice to select 'a' to 'b' along axis (numpy: [...,a:b,...])
        # else
        #     slice everything (numpy: [...,:,...])
        vix = tuple(slice(sms[i],ms[i]) if i in ssa else slice(None) for i in range(nd))

        # append values by copying
        mat[six] = self._mat    # copy self into matrix
        mat[vix] = values       # copy values into matrix
        self._mat = mat         # change pointer reference

        # append values
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

        if axis in self.square_axes:
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
        for axis in self.square_axes:
            self._mat = numpy.delete(self._mat, obj, axis = axis)
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
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    # TODO: FIXME
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
            raise ValueError("'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i not in self.square_axes) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))
        if (self._taxa is not None) and (taxa is None):
            taxa = numpy.empty(values.shape[self.taxa_axis], dtype = "object")  # fill with None
        if (self._taxa_grp is not None) and (taxa_grp is None):
            raise TypeError("cannot incorp: 'taxa_grp' argument is required")

        # TODO: # FIXME: figure out incorporation logic
        # incorp values
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
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
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
            raise ValueError(
                "cannot lexsort on axes {0} (taxa axes): {1}".format(
                    self.square_axes, emess
                )
            )

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
        Reorder the VariantMatrix.

        Parameters
        ----------
        indices : numpy.ndarray
            Indices of where to place elements.
        axis : int
            The axis over which to reorder values.
        """
        axis = get_axis(axis, self.mat_ndim)                   # transform axis number to an index

        if axis in self.square_axes:
            self.reorder_taxa(indices = indices, **kwargs)
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
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # reorder arrays
        for axis in self.square_axes:           # for each taxa axis
            ix = tuple(                         # build a tuple to slice the matrix
                indices if i == axis else slice(None) for i in range(self.mat_ndim)
            )
            self._mat = self._mat[ix]           # reorder matrix
        if self._taxa is not None:                      # if we have taxa names
            self._taxa = self._taxa[indices]            # reorder taxa array
        if self._taxa_grp is not None:                  # if we have taxa groups
            self._taxa_grp = self._taxa_grp[indices]    # reorder taxa group array

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
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis in self.square_axes:
            self.sort_taxa(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    # sort_taxa is unaltered

    ################### Grouping Methods ###################
    def group(
            self, 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Sort matrix along axis, then populate grouping indices for the axis.
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis in self.square_axes:
            self.group_taxa(**kwargs)
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
        if axis in self.square_axes:
            self.ungroup_taxa(**kwargs)
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
        else:
            raise ValueError("cannot test for grouping along axis {0}".format(axis))

        return grouped

    # is_grouped_taxa is unaltered



################################## Utilities ###################################
def check_is_DenseSquareTaxaMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSquareTaxaMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSquareTaxaMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseSquareTaxaMatrix.__name__,type(v).__name__))
