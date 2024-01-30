"""
Module implementing a dense trait matrix with axes that are square and 
associated error checking routines.
"""

__all__ = [
    "DenseSquareTraitMatrix",
    "check_is_DenseSquareTraitMatrix",
]

import copy
from typing import Optional, Sequence, Union
import numpy
from numpy.typing import ArrayLike
from pybrops.core.error.error_attr_python import check_is_iterable
from pybrops.core.error.error_generic_python import generic_check_isinstance
from pybrops.core.error.error_type_python import check_is_array_like
from pybrops.core.mat.DenseSquareMatrix import DenseSquareMatrix
from pybrops.core.mat.DenseTraitMatrix import DenseTraitMatrix
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.SquareTraitMatrix import SquareTraitMatrix
from pybrops.core.mat.util import get_axis

class DenseSquareTraitMatrix(
        DenseSquareMatrix,
        DenseTraitMatrix,
        SquareTraitMatrix,
    ):
    """
    A concrete class for dense matrices with trait metadata and axes that are
    square.

    The purpose of this abstract class is to merge the following implementations
    and interfaces:

        1. DenseSquareMatrix (implementation)
        2. DenseTraitMatrix (implementation)
        3. SquareTraitMatrix (interface)
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            mat: numpy.ndarray,
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseSquareTraitMatrix.
        
        Parameters
        ----------
        mat : numpy.ndarray
            A numpy.ndarray used to construct the object.
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseSquareTraitMatrix, self).__init__(
            mat, 
            **kwargs
        )
        self.trait = trait

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseSquareTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A shallow copy of the matrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            trait = copy.copy(self.trait),
        )

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseSquareTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A deep copy of the matrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            trait = copy.deepcopy(self.trait),
        )

        return out

    ############################ Object Properties #############################

    ##################### Matrix Data ######################
    # mat inherited from DenseSquareMatrix

    ############## Square Metadata Properties ##############
    @property
    def nsquare(self) -> int:
        """Number of axes that are square."""
        return self.nsquare_trait

    @property
    def square_axes(self) -> tuple:
        """Axis indices for axes that are square."""
        return self.square_trait_axes

    @property
    def square_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        return self.square_trait_axes_len

    @property
    def nsquare_trait(self) -> int:
        """Number of trait axes that are square."""
        return len(self.square_trait_axes)
    
    @property
    def square_trait_axes(self) -> tuple:
        """Axis indices for trait axes that are square."""
        return (0,1)
    
    @property
    def square_trait_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        return tuple(self.mat_shape[ix] for ix in self.square_trait_axes)

    ################ Trait Data Properites #################

    ############## Trait Metadata Properites ###############
    @property
    def trait_axis(self) -> int:
        """First square axis along which trait are stored"""
        return min(self.square_trait_axes)

    ################### Fill data lookup ###################
    # _fill_value inherited from DenseSquareMatrix

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
        return self.is_square_trait()

    def is_square_trait(
            self
        ) -> bool:
        """
        Determine whether the trait axes lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square trait axes are the same length.
            ``False`` if not all square trait axes are the same length.
        """
        siter = iter(self.square_trait_axes_len)    # get iterator of axis lengths
        try:                                        # try next
            e0 = next(siter)                        # get the first element of siter
        except StopIteration:                       # catch StopIteration exception
            return False                            # something is horribly wrong; no axes!
        out = all(e0 == e for e in siter)           # determine if all are equal to e0
        return out

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
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
        trait : numpy.ndarray
            Trait names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A copy of DenseSquareTraitMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_trait_axes:
            out = self.adjoin_trait(
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_trait(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the Matrix.
        trait : numpy.ndarray
            Trait names to adjoin to the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : DenseSquareTraitMatrix
            A copy of DenseSquareTraitMatrix with values appended to axis. Note that adjoin does
            not occur in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i not in self.square_trait_axes) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            trait = numpy.empty(values.shape[self.trait_axis], dtype = "object")   # fill with None

        # calculate shape of new matrix to allocate
        mshape = list(self.mat_shape)       # get matrix shape of self
        vshape = values.shape               # get matrix shape of values
        for ix in self.square_trait_axes:   # for each square trait axis index
            mshape[ix] += vshape[ix]        # add value axis length to self shape
        mshape = tuple(mshape)              # convert list to tuple

        # allocate memory and fill with default fill values
        mat = numpy.empty(mshape, dtype = self._mat.dtype)
        mat.fill(self._fill_value[mat.dtype])

        # calculate indexing tuples for copying matrix values into 'mat'
        nd = mat.ndim           # number of dimensions
        sms = self.mat_shape    # matrix shape for self
        ssa = self.square_trait_axes  # tuple of square trait axes
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
        if self._trait is not None:
            trait = numpy.append(self._trait, trait, axis = 0)

        # create output
        out = self.__class__(
            mat = mat,
            trait = trait,
            **kwargs
        )

        return out

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
        if axis in self.square_trait_axes:
            out = self.delete_trait(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
        """
        Delete sub-arrays along the trait axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A DenseSquareTraitMatrix with deleted elements. Note that concat does not occur
            in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        # get values
        mat = self._mat
        trait = self._trait

        # delete values
        for axis in self.square_trait_axes:
            mat = numpy.delete(mat, obj, axis = axis)
        if trait is not None:
            trait = numpy.delete(trait, obj, axis = 0)

        # create output
        out = self.__class__(
            mat = mat,
            trait = trait,
            **kwargs
        )

        return out

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
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
        trait : numpy.ndarray
            Trait names to insert into the Matrix.
            If values is a DenseHaplotypeMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A DenseSquareTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis in self.square_trait_axes:
            out = self.insert_trait(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    # TODO: FIXME
    def insert_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
        """
        Insert values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : Matrix, numpy.ndarray
            Values to insert into the matrix.
        trait : numpy.ndarray
            Trait names to insert into the Matrix.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTraitMatrix
            A DenseSquareTraitMatrix with values inserted. Note that insert does not occur
            in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot insert: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i not in self.square_trait_axes) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            trait = numpy.empty(values.shape[self.trait_axis], dtype = "object")   # fill with None

        # TODO: # FIXME: figure out insertion logic
        # insert values
        values = numpy.insert(self._mat, obj, values, axis = self.trait_axis)
        if self._trait is not None:
            trait = numpy.insert(self._trait, obj, trait, axis = 0)

        # create output
        out = self.__class__(
            mat = values,
            trait = trait,
            **kwargs
        )

        return out

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
        if axis in self.square_trait_axes:
            out = self.select_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_trait(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
        """
        Select certain values from the Matrix along the trait axis.

        Parameters
        ----------
        indices : array_like (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseSquareTraitMatrix
            The output DenseSquareTraitMatrix with values selected. Note that select does not
            occur in-place: a new DenseSquareTraitMatrix is allocated and filled.
        """
        # check for array_like
        check_is_array_like(indices, "indices")

        # get values
        mat = self._mat
        trait = self._trait

        # select values
        for axis in self.square_trait_axes:
            mat = numpy.take(mat, indices, axis = axis)
        if trait is not None:
            trait = numpy.take(trait, indices, axis = 0)

        # create output
        out = self.__class__(
            mat = mat,
            trait = trait,
            **kwargs
        )

        return out

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
        if axis in mats[0].square_trait_axes:
            out = cls.concat_trait(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    # TODO: FIXME
    @classmethod
    def concat_trait(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseSquareTraitMatrix':
        """
        Concatenate list of Matrix together along the trait axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : DenseSquareTraitMatrix
            The concatenated DenseSquareTraitMatrix. Note that concat does not occur in-place:
            a new DenseSquareTraitMatrix is allocated and filled.
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

        # test matrix compatibility (same axis length along non-trait axes)
        for i,v in enumerate(shape_t):                              # for each index,tuple in shape_t
            if (i != mats[0].trait_axis) and any(l != v[0] for l in v): # if not the trait axis AND axis lengths are different
                raise ValueError("cannot concat: matrix shapes do not all align along axis {0}".format(i))

        # create matrix lists
        mat_ls = [m.mat for m in mats]
        trait_ls = [m.trait for m in mats]

        # process/error check trait_ls
        if all(e is None for e in trait_ls):                             # if all elements are None
            trait_ls = None                                              # replace list with None
        else:                                                           # else at least one matrix does not have a trait array
            for i,v in enumerate(trait_ls):                              # for each index,element in trait_ls
                if v is None:                                           # if element is None
                    ntrait = shape_t[mats[0].trait_axis][i]               # get number of trait
                    trait_ls[i] = numpy.empty(ntrait, dtype = "object")   # replace with array of None

        # concatenate mat, trait items
        mat = numpy.concatenate(mat_ls, axis = mats[0].trait_axis)
        trait = None if trait_ls is None else numpy.concatenate(trait_ls, axis = 0)

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseHaplotypeMatrix
        # use first element as source of variant data
        out = cls(
            mat = mat,
            trait = trait,
            **kwargs
        )

        return out

    ######### Matrix element in-place-manipulation #########
    def append(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
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
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis in self.square_trait_axes:
            self.append_trait(
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    def append_trait(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the Matrix along the trait axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        trait : numpy.ndarray
            Trait names to append to the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot append: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i not in self.square_trait_axes) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            trait = numpy.empty(values.shape[self.trait_axis], dtype = "object") # fill with None

        # calculate shape of new matrix to allocate
        mshape = list(self.mat_shape)       # get matrix shape of self
        vshape = values.shape               # get matrix shape of values
        for ix in self.square_trait_axes:   # for each square axis index
            mshape[ix] += vshape[ix]        # add value axis length to self shape
        mshape = tuple(mshape)              # convert list to tuple

        # allocate memory and fill with default fill values
        mat = numpy.empty(mshape, dtype = self._mat.dtype)
        mat.fill(self._fill_value[mat.dtype])

        # calculate indexing tuples for copying matrix values into 'mat'
        nd = mat.ndim                   # number of dimensions
        sms = self.mat_shape            # matrix shape for self
        ssa = self.square_trait_axes    # tuple of square trait axes
        ms = mat.shape                  # output matrix shape

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
        if self._trait is not None:
            self._trait = numpy.append(self._trait, trait, axis = 0)

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

        if axis in self.square_trait_axes:
            self.remove_trait(obj = obj, **kwargs)
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def remove_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the trait axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # delete values
        for axis in self.square_trait_axes:
            self._mat = numpy.delete(self._mat, obj, axis = axis)
        if self._trait is not None:
            self._trait = numpy.delete(self._trait, obj, axis = 0)

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
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
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        if axis in self.square_trait_axes:
            self.incorp_trait(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    # TODO: FIXME
    def incorp_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Incorporate values along the trait axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : Matrix, numpy.ndarray
            Values to incorporate into the matrix.
        trait : numpy.ndarray
            Trait names to incorporate into the Matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseHaplotypeMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i not in self.square_trait_axes) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            trait = numpy.empty(values.shape[self.trait_axis], dtype = "object")  # fill with None

        # TODO: # FIXME: figure out incorporation logic
        # incorp values
        self._mat = numpy.insert(self._mat, obj, values, axis = self.trait_axis)

        if self._trait is not None:
            self._trait = numpy.insert(self._trait, obj, trait, axis = 0)

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
        if axis in self.square_trait_axes:
            indices = self.lexsort_trait(keys = keys, **kwargs)
        else:
            raise ValueError("cannot lexsort along axis {0}".format(axis))

        return indices

    def lexsort_trait(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys along the trait
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
            keys = (self._trait,)         # trait default keys
            emess = "trait is None"           # trait error message

        # remove None keys
        keys = tuple(k for k in keys if k is not None)

        # raise error if no keys remain
        if len(keys) == 0:
            raise ValueError("cannot lexsort on axis {0} (trait axis): {1}".format(self.trait_axis, emess))

        # raise error if keys are of incompatible length
        if any(len(k) != self.ntrait for k in keys):
            emess = "keys are not all length {0}".format(self.ntrait)
            raise ValueError(
                "cannot lexsort on axes {0} (trait axes): {1}".format(
                    self.square_trait_axes, emess
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

        if axis in self.square_trait_axes:
            self.reorder_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot reorder along axis {0}".format(axis))

    def reorder_trait(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Reorder elements of the Matrix along the trait axis using an array of
        indices. Note this modifies the Matrix in-place.

        Parameters
        ----------
        indices : (N,) ndarray of ints
            Array of indices that reorder the matrix along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # reorder arrays
        for axis in self.square_trait_axes:           # for each trait axis
            ix = tuple(                         # build a tuple to slice the matrix
                indices if i == axis else slice(None) for i in range(self.mat_ndim)
            )
            self._mat = self._mat[ix]           # reorder matrix
        if self._trait is not None:                      # if we have trait names
            self._trait = self._trait[indices]            # reorder trait array

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
        if axis in self.square_trait_axes:
            self.sort_trait(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    # sort_trait is unaltered

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_DenseSquareTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type DenseSquareTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseSquareTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,DenseSquareTraitMatrix.__name__,type(v).__name__))
