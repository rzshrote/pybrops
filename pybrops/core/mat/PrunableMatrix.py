"""
Module defining interfaces and associated error checking routines for matrices
that can be pruned along an axis.
"""

from typing import Any

import numpy
from pybrops.core.mat.Matrix import Matrix

# TODO: is this class even necessary?
class PrunableMatrix(Matrix):
    """
    An abstract class for prunable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix column and row pruning
        2) Provide backwards compatibility for previous software iterations.

    The shape of a PrunableMatrix should be immutable.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        PrunableMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(PrunableMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def prune(
            self, 
            axis: int, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate a set of indices for selection of rows, columns, etc.

        Parameters
        ----------
        axis : int
            The axis along which indices are selected.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of indices from which selections can be made.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PrunableMatrix(v: Any) -> bool:
    """
    Determine whether an object is a PrunableMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PrunableMatrix object instance.
    """
    return isinstance(v, PrunableMatrix)

def check_is_PrunableMatrix(v: Any, varname: str) -> None:
    """
    Check if object is of type PrunableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PrunableMatrix):
        raise TypeError("'%s' must be a PrunableMatrix." % varname)
