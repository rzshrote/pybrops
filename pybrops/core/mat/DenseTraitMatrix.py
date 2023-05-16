"""
Module implementing a dense matrix with trait metadata and associated error
checking routines.
"""

import numpy
import copy
from typing import Any, Sequence, Union
from typing import Optional
from numpy.typing import ArrayLike

from pybrops.core.error.error_attr_python import check_is_iterable
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_generic_python import generic_check_isinstance
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix

class DenseTraitMatrix(DenseMutableMatrix,TraitMatrix):
    """
    A concrete class for dense matrices with trait metadata.

    The purpose of this concrete class is to implement base functionality for:
        1) Dense matrix trait metadata.
        2) Dense matrix trait routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DenseTraitMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            A numpy.ndarray used to construct the object.
        trait : numpy.ndarray, None
            A numpy.ndarray of shape ``(t,)`` containing trait names.
            If ``None``, do not store any trait name information.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseTraitMatrix, self).__init__(
            mat = mat,
            **kwargs
        )
        self.trait = trait

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DenseTraitMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DenseTraitMatrix
            A copy of the DenseTraitMatrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            trait = copy.copy(self.trait)
        )

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DenseTraitMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Deep copy metadata.

        Returns
        -------
        out : Matrix
            A deep copy of the DenseTraitMatrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            trait = copy.deepcopy(self.trait, memo)
        )

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ###################### Trait data ######################
    @property
    def trait(self) -> Union[numpy.ndarray,None]:
        """Trait label."""
        return self._trait
    @trait.setter
    def trait(self, value: Union[numpy.ndarray,None]) -> None:
        """Set trait label array"""
        if value is not None:
            check_is_ndarray(value, "trait")
            check_ndarray_dtype_is_object(value, "trait")
            check_ndarray_ndim(value, "trait", 1)
            check_ndarray_axis_len(value, "trait", 0, self.ntrait)
        self._trait = value
    @trait.deleter
    def trait(self) -> None:
        """Delete trait label array"""
        del self._trait
    
    #################### Trait metadata ####################
    @property
    def ntrait(self) -> int:
        """Number of traits."""
        return self._mat.shape[self.trait_axis]
    @ntrait.setter
    def ntrait(self, value: int) -> None:
        """Set number of traits"""
        error_readonly("ntrait")
    @ntrait.deleter
    def ntrait(self) -> None:
        """Delete number of traits"""
        error_readonly("ntrait")
    
    @property
    def trait_axis(self) -> int:
        """Axis along which traits are stored."""
        return 0
    @trait_axis.setter
    def trait_axis(self, value: int) -> None:
        """Set trait axis number"""
        error_readonly("ntrait")
    @trait_axis.deleter
    def trait_axis(self) -> None:
        """Delete trait axis number"""
        error_readonly("ntrait")
    
    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTraitMatrix':
        """
        Add additional elements to the end of the Matrix along an axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values to be adjoined to the Matrix.
        axis : int
            The axis along which values are adjoined.
        trait : numpy.ndarray
            Trait names to adjoin to the Matrix.
            If values is a DenseTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTraitMatrix
            A copy of mat with values adjoined to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.trait_axis:
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
        ) -> 'DenseTraitMatrix':
        """
        Add additional elements to the end of the Matrix along the trait axis.

        Parameters
        ----------
        values : DenseTraitMatrix, numpy.ndarray
            Values to be adjoined to the Matrix.
        trait : numpy.ndarray
            Trait names to adjoin to the Matrix.
            If values is a DenseTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTraitMatrix
            A copy of mat with values adjoined to axis. Note that adjoin does
            not occur in-place: a new Matrix is allocated and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__): # TODO: change to DenseTraitMatrix instead of self.__class__???
            if trait is None:
                trait = values.trait
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.trait_axis) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            raise TypeError("cannot adjoin: 'trait' argument is required")

        # adjoin values
        values = numpy.append(self._mat, values, axis = self.trait_axis)
        if self._trait is not None:
            trait = numpy.append(self._trait, trait, axis = 0)

        out = self.__class__(
            mat = values,
            trait = trait,
            **kwargs
        )

        return out

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.trait_axis:
            out = self.delete_trait(
                obj = obj,
                **kwargs
            )
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get values
        mat = self._mat
        trait = self._trait

        # delete values
        mat = numpy.delete(mat, obj, axis = self.trait_axis)
        if trait is not None:
            trait = numpy.delete(trait, obj, axis = 0)

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
        ) -> 'DenseTraitMatrix':
        """
        Insert values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            inserted.
        values : DenseTraitMatrix, numpy.ndarray
            Values to insert into the matrix.
        axis : int
            The axis along which values are inserted.
        trait : numpy.ndarray
            Trait names to insert into the Matrix.
            If values is a DenseTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseTraitMatrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.trait_axis:
            out = self.insert_trait(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    def insert_trait(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            trait: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            A Matrix with values inserted. Note that insert does not occur
            in-place: a new Matrix is allocated and filled.
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
            if (i != self.trait_axis) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            raise ValueError("cannot insert: 'trait' argument is required")

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
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            The output matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.trait_axis:
            out = self.select_trait(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_trait(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            The output Matrix with values selected. Note that select does not
            occur in-place: a new Matrix is allocated and filled.
        """
        # get values
        mat = self._mat
        trait = self._trait

        # select values
        mat = numpy.take(mat, indices, axis = self.trait_axis)
        if trait is not None:
            trait = numpy.take(trait, indices, axis = 0)

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
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].trait_axis:
            out = cls.concat_trait(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_trait(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DenseTraitMatrix':
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
        out : DenseTraitMatrix
            The concatenated DenseTraitMatrix. Note that concat does not occur in-place:
            a new DenseTraitMatrix is allocated and filled.
        """
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # ensure that we have an array_like of length >= 1
        if len(mats) <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DenseTraitMatrix
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
        if all(e is None for e in trait_ls):    # if all elements are None
            trait_ls = None                     # replace list with None
        elif any(e is None for e in trait_ls):  # else if any elements are None
            raise ValueError("cannot concat: 'trait' needed for all Matrix in list")

        # concatenate mat, trait items
        mat = numpy.concatenate(mat_ls, axis = mats[0].trait_axis)
        trait = None if trait_ls is None else numpy.concatenate(trait_ls, axis = 0)

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseTraitMatrix
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
        values : DenseTraitMatrix, numpy.ndarray
            Values are appended to append to the matrix.
            Must be of type int8.
            Must be of shape (m, n, p)
        axis : int
            The axis along which values are appended.
        trait : numpy.ndarray
            Trait names to append to the Matrix.
            If values is a DenseTraitMatrix that has a non-None
            trait field, providing this argument overwrites the field.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.trait_axis:
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
            if (i != self.trait_axis) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            raise ValueError("cannot append: 'trait' argument is required")

        # append values
        self._mat = numpy.append(self._mat, values, axis = self.trait_axis)
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

        if axis == self.trait_axis:
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
        self._mat = numpy.delete(self._mat, obj, axis = self.trait_axis)

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

        if axis == self.trait_axis:
            self.incorp(
                obj = obj,
                values = values,
                trait = trait,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

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
            raise ValueError("'values' must be of type DenseTraitMatrix or numpy.ndarray")

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.trait_axis) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))
        if (self._trait is not None) and (trait is None):
            raise ValueError("cannot incorp: 'trait' argument is required")

        # insert values
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
        if axis == self.trait_axis:
            self.lexsort_trait(keys = keys, **kwargs)
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
            keys = (self._trait,)    # trait default keys
            emess = "trait is None"  # trait error message

        # remove None keys
        keys = tuple(k for k in keys if k is not None)

        # raise error if no keys remain
        if len(keys) == 0:
            raise ValueError("cannot lexsort on axis {0} (trait axis): {1}".format(self.trait_axis, emess))

        # raise error if keys are of incompatible length
        if any(len(k) != self.ntrait for k in keys):
            emess = "keys are not all length {0}".format(self.ntrait)
            raise ValueError("cannot lexsort on axis {0} (trait axis): {1}".format(self.trait_axis, emess))

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

        if axis == self.trait_axis:
            self.reorder(indices = indices, **kwargs)
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
        # build a tuple to slice the matrix
        ix = tuple(indices if i == self.trait_axis else slice(None) for i in range(self.mat_ndim))

        # reorder arrays
        self._mat = self._mat[ix]

        if self._trait is not None:
            self._trait = self._trait[indices]                # reorder trait array

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
        """
        # transform axis number to an index
        axis = get_axis(axis, self.mat_ndim)

        # dispatch functions
        if axis == self.trait_axis:
            self.sort_trait(keys = keys, **kwargs)
        else:
            raise ValueError("cannot sort along axis {0}".format(axis))

    def sort_trait(
            self, 
            keys: Union[tuple,numpy.ndarray,None] = None, 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the Matrix along the trait axis using a sequence of
        keys. Note this modifies the Matrix in-place.

        Parameters
        ----------
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        kwargs : dict
            Additional keyword arguments.
        """
        # get indices for sort
        indices = self.lexsort_trait(keys, **kwargs)

        # reorder internals
        self.reorder_trait(indices, **kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseTraitMatrix(v: Any) -> bool:
    """
    Determine whether an object is a DenseTraitMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DenseTraitMatrix object instance.
    """
    return isinstance(v, DenseTraitMatrix)

def check_is_DenseTraitMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type DenseTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_DenseTraitMatrix(v):
        raise TypeError("'{0}' must be a DenseTraitMatrix".format(vname))
