"""
Module defining interfaces and associated error checking routines for matrices
that can be pruned along an axis.
"""

__all__ = [
    "PrunableMatrix",
    "check_is_PrunableMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
import numpy
from pybrops.core.mat.Matrix import Matrix

# TODO: is this class even necessary?
class PrunableMatrix(Matrix,metaclass=ABCMeta):
    """
    An abstract class for prunable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix column and row pruning
        2) Provide backwards compatibility for previous software iterations.

    The shape of a PrunableMatrix should be immutable.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
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



################################## Utilities ###################################
def check_is_PrunableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type PrunableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PrunableMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,PrunableMatrix.__name__,type(v).__name__))
