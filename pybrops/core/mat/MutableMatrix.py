"""
Module defining mutable matrix interfaces and associated error checking routines.
"""

__all__ = [
    "MutableMatrix",
    "check_is_MutableMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Sequence
from typing import Union

import numpy
from pybrops.core.mat.Matrix import Matrix

class MutableMatrix(
        Matrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for mutable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix shape changing routines.

    The shape of a MutableMatrix is mutable.
    """

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__         inherited from ``Matrix``
    ### __sub__         inherited from ``Matrix``
    ### __mul__         inherited from ``Matrix``
    ### __matmul__      inherited from ``Matrix``
    ### __truediv__     inherited from ``Matrix``
    ### __floordiv__    inherited from ``Matrix``
    ### __mod__         inherited from ``Matrix``
    ### __divmod__      inherited from ``Matrix``
    ### __pow__         inherited from ``Matrix``
    ### __lshift__      inherited from ``Matrix``
    ### __rshift__      inherited from ``Matrix``
    ### __and__         inherited from ``Matrix``
    ### __xor__         inherited from ``Matrix``
    ### __or__          inherited from ``Matrix``

    ############# Backwards numeric operators ##############
    ### __radd__        inherited from ``Matrix``
    ### __rsub__        inherited from ``Matrix``
    ### __rmul__        inherited from ``Matrix``
    ### __rmatmul__     inherited from ``Matrix``
    ### __rtruediv__    inherited from ``Matrix``
    ### __rfloordiv__   inherited from ``Matrix``
    ### __rmod__        inherited from ``Matrix``
    ### __rdivmod__     inherited from ``Matrix``
    ### __rlshift__     inherited from ``Matrix``
    ### __rrshift__     inherited from ``Matrix``
    ### __rand__        inherited from ``Matrix``
    ### __rxor__        inherited from ``Matrix``
    ### __ror__         inherited from ``Matrix``

    ############# Augmented numeric operators ##############
    ### __iadd__        inherited from ``Matrix``
    ### __isub__        inherited from ``Matrix``
    ### __imul__        inherited from ``Matrix``
    ### __imatmul__     inherited from ``Matrix``
    ### __itruediv__    inherited from ``Matrix``
    ### __ifloordiv__   inherited from ``Matrix``
    ### __imod__        inherited from ``Matrix``
    ### __ipow__        inherited from ``Matrix``
    ### __ilshift__     inherited from ``Matrix``
    ### __irshift__     inherited from ``Matrix``
    ### __iand__        inherited from ``Matrix``
    ### __ixor__        inherited from ``Matrix``
    ### __ior__         inherited from ``Matrix``

    ################## Logical operators ###################
    ### __lt__          inherited from ``Matrix``
    ### __le__          inherited from ``Matrix``
    ### __eq__          inherited from ``Matrix``
    ### __ne__          inherited from ``Matrix``
    ### __gt__          inherited from ``Matrix``
    ### __ge__          inherited from ``Matrix``

    ################# Container operators ##################
    ### __len__         inherited from ``Matrix``
    ### __getitem__     inherited from ``Matrix``
    ### __setitem__     inherited from ``Matrix``
    ### __delitem__     inherited from ``Matrix``
    ### __iter__        inherited from ``Matrix``

    #################### Matrix copying ####################
    ### __copy__        inherited from ``Matrix``
    ### __deepcopy__    inherited from ``Matrix``

    ########### Miscellaneous special functions ############
    ### __repr__        inherited from ``Matrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat             inherited from ``Matrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim        inherited from ``Matrix``
    ### mat_shape       inherited from ``Matrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy            inherited from ``Matrix``
    ### deepcopy        inherited from ``Matrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin          inherited from ``Matrix``
    ### delete          inherited from ``Matrix``
    ### insert          inherited from ``Matrix``
    ### select          inherited from ``Matrix``

    ######### Matrix element in-place-manipulation #########
    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
    def incorp(
            self, 
            obj: Union[int,slice,Sequence], 
            values: Union[Matrix,numpy.ndarray], 
            axis: int, 
            **kwargs: dict
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

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat          inherited from ``Matrix``



################################## Utilities ###################################
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
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,MutableMatrix.__name__,type(v).__name__))
