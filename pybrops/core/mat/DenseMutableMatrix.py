"""
Module implementing a dense mutable matrix and associated error checking routines.

Mutable refers to the ability of adding/removing rows and columns from the matrix.
"""

from typing import Any, Sequence, Union
import numpy
from pybrops.core.mat.util import get_axis
from pybrops.core.mat.DenseMatrix import DenseMatrix
from pybrops.core.mat.MutableMatrix import MutableMatrix

class DenseMutableMatrix(DenseMatrix,MutableMatrix):
    """
    A concrete class for dense mutable matrices.
    Dense mutable matrices utilize numpy.ndarray's for data storage.

    The purpose of this concrete class is to implement base functionality for:
        1) Dense matrix in-place matrix shape changing routines.
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
        Constructor for DenseMutableMatrix

        Parameters
        ----------
        mat : numpy.ndarray
            Dense matrix to use to create object.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseMutableMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ######### Matrix element in-place-manipulation #########
    def append(
            self, 
            values: Union[DenseMatrix,numpy.ndarray], 
            axis: int = -1, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the matrix.

        Parameters
        ----------
        values : DenseMatrix, numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a DenseMatrix extract Matrix.mat values
        if isinstance(values, DenseMatrix):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # append values
        self._mat = numpy.append(self._mat, values, axis)

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
        obj : slice, int, or Sequence of ints
            Indicate indices of sub-arrays to remove along the specified axis.
        axis: int
            The axis along which to remove the subarray defined by obj.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # delete values
        self._mat = numpy.delete(self._mat, obj, axis)

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[DenseMatrix,numpy.ndarray], 
            axis: int = -1, 
            **kwargs
        ) -> None:
        """
        Incorporate values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : DenseMatrix or numpy.ndarray
            Values to incorporate into the matrix.
        axis : int
            The axis along which values are incorporated.
        kwargs : dict
            Additional keyword arguments.
        """
        # get axis
        axis = get_axis(axis, self._mat.ndim)

        # if given a Matrix extract Matrix.mat values
        if isinstance(values, DenseMatrix):
            values = values.mat
        elif not isinstance(values, numpy.ndarray):
            raise ValueError("'values' must be of type DenseMatrix or numpy.ndarray")

        # incorporate values
        self._mat = numpy.insert(self._mat, obj, values, axis)



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseMutableMatrix(v: object) -> bool:
    """
    Determine whether an object is a HDF5InputOutput.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a HDF5InputOutput object instance.
    """
    return isinstance(v, DenseMutableMatrix)

def check_is_DenseMutableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type HDF5InputOutput. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseMutableMatrix):
        raise TypeError("'%s' must be a DenseMutableMatrix." % varname)
