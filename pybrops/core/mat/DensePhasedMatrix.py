"""
Module implementing a dense phased matrix and associated error checking routines.
"""

import copy
import numpy
from typing import Any, Optional, Sequence, Union
from numpy.typing import ArrayLike

from pybrops.core.error import check_is_array_like
from pybrops.core.error import check_is_iterable
from pybrops.core.error import error_readonly
from pybrops.core.error import generic_check_isinstance
from pybrops.core.mat.Matrix import Matrix
from pybrops.core.mat.util import get_axis
from pybrops.core.mat.DenseMutableMatrix import DenseMutableMatrix
from pybrops.core.mat.PhasedMatrix import PhasedMatrix

class DensePhasedMatrix(DenseMutableMatrix,PhasedMatrix):
    """
    A concrete class implementing dense phased matrices.
    A phased matrix is defined as a matrix with a third dimension.

    Dense phased matrices utilize numpy.ndarray's for data storage.

    The purpose of this concrete class is to implement base functionality for:
        1) Dense matrix phase manipulation routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the concrete class DensePhasedMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Matrix used to construct the object.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DensePhasedMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DensePhasedMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
        )

        return out

    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DensePhasedMatrix':
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
        )

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Phase Metadata Properites ###############
    def nphase():
        doc = "Number of chromosome phases represented by the matrix."
        def fget(self):
            """Get number of phases"""
            return self._mat.shape[self.phase_axis]
        def fset(self, value):
            """Set number of phases"""
            error_readonly("nphase")
        def fdel(self):
            """Delete number of phases"""
            error_readonly("nphase")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    nphase = property(**nphase())

    def phase_axis():
        doc = "Axis along which phases are stored property."
        def fget(self):
            """Get phase axis number"""
            return 0
        def fset(self, value):
            """Set phase axis number"""
            error_readonly("phase_axis")
        def fdel(self):
            """Delete phase axis number"""
            error_readonly("phase_axis")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    phase_axis = property(**phase_axis())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element copy-on-manipulation ##########
    def adjoin(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
        """
        Add additional elements to the end of the DensePhasedMatrix along an axis.

        Parameters
        ----------
        values : DensePhasedMatrix, numpy.ndarray
            Values are appended to append to the Matrix.
        axis : int
            The axis along which values are adjoined.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedMatrix
            A copy of mat with values appended to axis. Note that adjoin does
            not occur in-place: a new DensePhasedMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.adjoin_phase(
                values = values,
                **kwargs
            )
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

        return out

    def adjoin_phase(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
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
        out : DensePhasedMatrix
            A copy of the DensePhasedMatrix with values adjoined along the phase axis.
            Note that adjoin does not occur in-place: a new DensePhasedMatrix is allocated and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("cannot adjoin: 'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot adjoin: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot adjoin: axis lengths incompatible for axis {0}".format(i))

        # adjoin values
        values = numpy.append(self._mat, values, axis = self.phase_axis)

        out = self.__class__(
            mat = values,
            **kwargs
        )

        return out

    def delete(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
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
        out : DensePhasedMatrix
            A DensePhasedMatrix with deleted elements. Note that concat does not occur
            in-place: a new DensePhasedMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.delete_phase(obj = obj, **kwargs)
        else:
            raise ValueError("cannot delete along axis {0}".format(axis))

        return out

    def delete_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
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
        out : Matrix
            A Matrix with deleted elements. Note that concat does not occur
            in-place: a new Matrix is allocated and filled.
        """
        # get values
        mat = self._mat

        # delete values
        mat = numpy.delete(mat, obj, axis = self.phase_axis)

        out = self.__class__(
            mat = mat,
            **kwargs
        )

        return out

    def insert(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedMatrix
            A DensePhasedMatrix with values inserted. Note that insert does not occur
            in-place: a new DensePhasedMatrix is allocated and filled.
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
        else:
            raise ValueError("cannot insert along axis {0}".format(axis))

        return out

    def insert_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
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
        out : DensePhasedMatrix
            A DensePhasedMatrix with values inserted. Note that insert does not occur
            in-place: a new DensePhasedMatrix is allocated and filled.
        """
        # extract mat values
        if isinstance(values, self.__class__):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot insert: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot insert: axis lengths incompatible for axis {0}".format(i))

        # insert values
        values = numpy.insert(self._mat, obj, values, axis = self.phase_axis)

        # create output
        out = self.__class__(
            mat = values,
            **kwargs
        )

        return out

    def select(
            self, 
            indices: ArrayLike, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
        """
        Select certain values from the DensePhasedMatrix.

        Parameters
        ----------
        indices : ArrayLike (Nj, ...)
            The indices of the values to select.
        axis : int
            The axis along which values are selected.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedMatrix
            The output DensePhasedMatrix with values selected. Note that select does not
            occur in-place: a new DensePhasedMatrix is allocated and filled.
        """
        axis = get_axis(axis, self.mat_ndim)    # get axis
        out = None                              # declare variable

        # dispatch functions to handle operations
        if axis == self.phase_axis:
            out = self.select_phase(indices = indices, **kwargs)
        else:
            raise ValueError("cannot select along axis {0}".format(axis))

        return out

    def select_phase(
            self, 
            indices: ArrayLike, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
        """
        Select certain values from the DensePhasedMatrix along the phase axis.

        Parameters
        ----------
        indices : ArrayLike (Nj, ...)
            The indices of the values to select.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DensePhasedMatrix
            The output DensePhasedMatrix with values selected. Note that select does not
            occur in-place: a new DensePhasedMatrix is allocated and filled.
        """
        # check for array_like
        check_is_array_like(indices, "indices")

        # select values
        mat = numpy.take(self._mat, indices, axis = self.phase_axis)

        out = self.__class__(
            mat = mat,
            **kwargs
        )

        return out

    @classmethod
    def concat(
            cls, 
            mats: Sequence, 
            axis: int = -1, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
        """
        Concatenate a sequence of Matrix together along an axis.

        Parameters
        ----------
        mats : Sequence of Matrix
            List of Matrix to concatenate. The matrices must have the same
            shape, except in the dimension corresponding to axis.
        axis : int
            The axis along which the arrays will be joined.
        kwargs : dict
            Additional keyword arguments

        Returns
        -------
        out : DensePhasedMatrix
            The concatenated DensePhasedMatrix. Note that concat does not occur in-place:
            a new DensePhasedMatrix is allocated and filled.
        """
        axis = get_axis(axis, mats[0].mat_ndim)     # get axis
        out = None                                  # declare variable

        # dispatch items to worker functions
        if axis == mats[0].phase_axis:
            out = cls.concat_phase(mats, **kwargs)
        else:
            raise ValueError("cannot concat along axis {0}".format(axis))

        return out

    @classmethod
    def concat_phase(
            cls, 
            mats: Sequence, 
            **kwargs: dict
        ) -> 'DensePhasedMatrix':
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
        out : Matrix
            The concatenated matrix. Note that concat does not occur in-place:
            a new Matrix is allocated and filled.
        """
        # ensure that we have an iterable object
        check_is_iterable(mats, "mats")

        # ensure that we have an array_like of length >= 1
        if len(mats) <= 0:
            raise ValueError("need at least one Matrix to concatenate")

        # ensure that all items in mats are DensePhasedMatrix
        for i,v in enumerate(mats):
            generic_check_isinstance(v, "mats[{0}]".format(i), cls)

        # make sure dimensions are all identical to first element in mats
        if any(m.mat_ndim != mats[0].mat_ndim for m in mats):
            raise ValueError("cannot concat: not all matrices have the same number of dimensions")

        # extract tuple of shapes for testing compatibility
        shape_t = tuple(zip(*[m.mat.shape for m in mats]))

        # test matrix compatibility (same axis length along non-taxa axes)
        for i,v in enumerate(shape_t):                              # for each index,tuple in shape_t
            if (i != mats[0].phase_axis) and any(l != v[0] for l in v): # if not the taxa axis AND axis lengths are different
                raise ValueError("cannot concat: matrix shapes do not all align along axis {0}".format(i))

        # create matrix lists
        mat_ls = [m.mat for m in mats]

        # concatenate items
        mat = numpy.concatenate(mat_ls, axis = mats[0].phase_axis)

        # TODO: decide if first element in list is good source of information
        # concatenate everything and put into new DenseHaplotypeMatrix
        # use first element as source of variant data
        out = cls(
            mat = mat,
            **kwargs
        )

        return out

    ######### Matrix element in-place-manipulation #########
    def append(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int = -1, 
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
        if axis == self.phase_axis:
            self.append_phase(values, **kwargs)
        else:
            raise ValueError("cannot append along axis {0}".format(axis))

    def append_phase(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Append values to the Matrix along the phase axis.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        kwargs : dict
            Additional keyword arguments.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if isinstance(values, self.__class__):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot append: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot append: axis lengths incompatible for axis {0}".format(i))

        # append values
        self._mat = numpy.append(self._mat, values, axis = self.phase_axis)

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
        else:
            raise ValueError("cannot remove along axis {0}".format(axis))

    def remove_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            **kwargs: dict
        ) -> None:
        """
        Remove sub-arrays along the phase axis.

        Parameters
        ----------
        obj : int, slice, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        kwargs : dict
            Additional keyword arguments.
        """
        # delete values
        self._mat = numpy.delete(self._mat, obj, axis = self.phase_axis)

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: ArrayLike, 
            axis: int = -1, 
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
            self.incorp(
                obj = obj,
                values = values,
                **kwargs
            )
        else:
            raise ValueError("cannot incorp along axis {0}".format(axis))

    def incorp_phase(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
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
        kwargs : dict
            Additional keyword arguments.
        """
        # if given a DensePhasedGenotypeMatrix extract *.mat values
        if isinstance(values, self.__class__):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type {0} or numpy.ndarray".format(self.__class__))

        # perform error checks before allocating memory
        if values.ndim != self.mat_ndim:
            raise ValueError("cannot incorp: 'values' must have ndim == {0}".format(self.mat_ndim))
        for i,(j,k) in enumerate(zip(values.shape, self.mat_shape)):
            if (i != self.phase_axis) and (j != k):
                raise ValueError("cannot incorp: axis lengths incompatible for axis {0}".format(i))

        # insert values
        self._mat = numpy.insert(self._mat, obj, values, axis = self.phase_axis)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedMatrix(v: Any) -> bool:
    """
    Determine whether an object is a DensePhasedMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a DensePhasedMatrix object instance.
    """
    return isinstance(v, DensePhasedMatrix)

def check_is_DensePhasedMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type DensePhasedMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DensePhasedMatrix):
        raise TypeError("'%s' must be a DensePhasedMatrix." % vname)
