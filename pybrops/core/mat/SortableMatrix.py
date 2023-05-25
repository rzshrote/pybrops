"""
Module defining interfaces and associated error checking routines for matrices
that can have their axes sorted.
"""

from typing import Any, Sequence, Union

import numpy
from pybrops.core.mat.MutableMatrix import MutableMatrix

class SortableMatrix(MutableMatrix):
    """
    An abstract class for sortable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix sorting routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        SortableMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(SortableMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
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
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        axis : int
            Axis to be indirectly sorted.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        indices : (N,) ndarray of ints
            Array of indices that sort the keys along the specified axis.
        """
        raise NotImplementedError("method is abstract")

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
        indices : (N,) ndarray of ints, Sequence of ints
            Array of indices that reorder the matrix along the specified axis.
        axis : int
            Axis to be reordered.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

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
        keys : (k, N) array or tuple containing k (N,)-shaped sequences
            The k different columns to be sorted. The last column (or row if
            keys is a 2D array) is the primary sort key.
        axis : int
            Axis to be indirectly sorted.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_SortableMatrix(v: object) -> bool:
    """
    Determine whether an object is a SortableMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a SortableMatrix object instance.
    """
    return isinstance(v, SortableMatrix)

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
        raise TypeError("'%s' must be a SortableMatrix." % vname)
