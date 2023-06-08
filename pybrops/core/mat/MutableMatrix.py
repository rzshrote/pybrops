"""
Module defining mutable matrix interfaces and associated error checking routines.
"""

from typing import Any, Sequence, Union

import numpy
from pybrops.core.mat.Matrix import Matrix

class MutableMatrix(Matrix):
    """
    An abstract class for mutable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix shape changing routines.

    The shape of a MutableMatrix is mutable.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        MutableMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(MutableMatrix, self).__init__(**kwargs)

    ############################## Object Methods ##############################

    ######### Matrix element in-place-manipulation #########
    def append(
            self, 
            values: Union[Matrix,numpy.ndarray], 
            axis: int, 
            **kwargs: dict
        ) -> None:
        """
        Append values to the Matrix.

        Parameters
        ----------
        values : Matrix, numpy.ndarray
            Values are appended to append to the matrix.
        axis : int
            The axis along which values are appended.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def remove(
            self, 
            obj: Union[int,slice,Sequence], 
            axis: int, 
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
        raise NotImplementedError("method is abstract")

    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int, 
            **kwargs
        ) -> None:
        """
        Incorporate values along the given axis before the given indices.

        Parameters
        ----------
        obj: int, slice, or Sequence of ints
            Object that defines the index or indices before which values is
            incorporated.
        values : numpy.ndarray
            Values to incorporate into the matrix.
        axis : int
            The axis along which values are incorporated.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_MutableMatrix(v: object) -> bool:
    """
    Determine whether an object is a MutableMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a MutableMatrix object instance.
    """
    return isinstance(v, MutableMatrix)

def check_is_MutableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type MutableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MutableMatrix):
        raise TypeError("'%s' must be a MutableMatrix." % vname)
