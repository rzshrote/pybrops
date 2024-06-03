"""
Module defining interfaces and associated error checking routines for matrices
that can have their axes sorted.
"""

__all__ = [
    "SortableMatrix",
    "check_is_SortableMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from typing import Sequence
from typing import Union

import numpy
from pybrops.core.mat.MutableMatrix import MutableMatrix

class SortableMatrix(
        MutableMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for sortable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix sorting routines.
    """

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__         inherited from ``MutableMatrix``
    ### __sub__         inherited from ``MutableMatrix``
    ### __mul__         inherited from ``MutableMatrix``
    ### __matmul__      inherited from ``MutableMatrix``
    ### __truediv__     inherited from ``MutableMatrix``
    ### __floordiv__    inherited from ``MutableMatrix``
    ### __mod__         inherited from ``MutableMatrix``
    ### __divmod__      inherited from ``MutableMatrix``
    ### __pow__         inherited from ``MutableMatrix``
    ### __lshift__      inherited from ``MutableMatrix``
    ### __rshift__      inherited from ``MutableMatrix``
    ### __and__         inherited from ``MutableMatrix``
    ### __xor__         inherited from ``MutableMatrix``
    ### __or__          inherited from ``MutableMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__        inherited from ``MutableMatrix``
    ### __rsub__        inherited from ``MutableMatrix``
    ### __rmul__        inherited from ``MutableMatrix``
    ### __rmatmul__     inherited from ``MutableMatrix``
    ### __rtruediv__    inherited from ``MutableMatrix``
    ### __rfloordiv__   inherited from ``MutableMatrix``
    ### __rmod__        inherited from ``MutableMatrix``
    ### __rdivmod__     inherited from ``MutableMatrix``
    ### __rlshift__     inherited from ``MutableMatrix``
    ### __rrshift__     inherited from ``MutableMatrix``
    ### __rand__        inherited from ``MutableMatrix``
    ### __rxor__        inherited from ``MutableMatrix``
    ### __ror__         inherited from ``MutableMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__        inherited from ``MutableMatrix``
    ### __isub__        inherited from ``MutableMatrix``
    ### __imul__        inherited from ``MutableMatrix``
    ### __imatmul__     inherited from ``MutableMatrix``
    ### __itruediv__    inherited from ``MutableMatrix``
    ### __ifloordiv__   inherited from ``MutableMatrix``
    ### __imod__        inherited from ``MutableMatrix``
    ### __ipow__        inherited from ``MutableMatrix``
    ### __ilshift__     inherited from ``MutableMatrix``
    ### __irshift__     inherited from ``MutableMatrix``
    ### __iand__        inherited from ``MutableMatrix``
    ### __ixor__        inherited from ``MutableMatrix``
    ### __ior__         inherited from ``MutableMatrix``

    ################## Logical operators ###################
    ### __lt__          inherited from ``MutableMatrix``
    ### __le__          inherited from ``MutableMatrix``
    ### __eq__          inherited from ``MutableMatrix``
    ### __ne__          inherited from ``MutableMatrix``
    ### __gt__          inherited from ``MutableMatrix``
    ### __ge__          inherited from ``MutableMatrix``

    ################# Container operators ##################
    ### __len__         inherited from ``MutableMatrix``
    ### __getitem__     inherited from ``MutableMatrix``
    ### __setitem__     inherited from ``MutableMatrix``
    ### __delitem__     inherited from ``MutableMatrix``
    ### __iter__        inherited from ``MutableMatrix``

    #################### Matrix copying ####################
    ### __copy__        inherited from ``MutableMatrix``
    ### __deepcopy__    inherited from ``MutableMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__        inherited from ``MutableMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat             inherited from ``MutableMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim        inherited from ``MutableMatrix``
    ### mat_shape       inherited from ``MutableMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy            inherited from ``MutableMatrix``
    ### deepcopy        inherited from ``MutableMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin          inherited from ``MutableMatrix``
    ### delete          inherited from ``MutableMatrix``
    ### insert          inherited from ``MutableMatrix``
    ### select          inherited from ``MutableMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append          inherited from ``MutableMatrix``
    ### remove          inherited from ``MutableMatrix``
    ### incorp          inherited from ``MutableMatrix``

    ################### Sorting Methods ####################
    @abstractmethod
    def lexsort(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            axis: int, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Perform an indirect stable sort using a sequence of keys.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        axis : int
            Axis to be indirectly sorted.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : A (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def reorder(
            self, 
            indices: Union[numpy.ndarray,Sequence], 
            axis: int, 
            **kwargs: dict
        ) -> None:
        """
        Reorder elements of the Matrix using an array of indices. Note this
        modifies the Matrix in-place.

        Parameters
        ----------
        indices : A (N,) ndarray of ints, Sequence of ints
            Array of indices that reorder the matrix along the specified axis.
        axis : int
            Axis to be reordered.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def sort(
            self, 
            keys: Union[tuple,numpy.ndarray], 
            axis: int, 
            **kwargs: dict
        ) -> None:
        """
        Sort slements of the Matrix using a sequence of keys. Note this modifies
        the Matrix in-place.

        Parameters
        ----------
        keys : A (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        axis : int
            Axis to be indirectly sorted.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat          inherited from ``MutableMatrix``



################################## Utilities ###################################
def check_is_SortableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type SortableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SortableMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SortableMatrix.__name__,type(v).__name__))
