"""
Module defining interfaces and associated error checking routines for matrices
that can have their axes grouped.
"""

__all__ = [
    "GroupableMatrix",
    "check_is_GroupableMatrix",
]

from abc import ABCMeta, abstractmethod
from pybrops.core.mat.SortableMatrix import SortableMatrix

class GroupableMatrix(SortableMatrix,metaclass=ABCMeta):
    """
    An abstract class for groupable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix axis grouping routines.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################

    ################### Grouping Methods ###################
    @abstractmethod
    def group(
            self, 
            axis: int, 
            **kwargs: dict
        ) -> None:
        """
        Sort the GroupableMatrix along an axis, then populate grouping indices.

        Parameters
        ----------
        axis : int
            The axis along which values are grouped.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def ungroup(
            self,
            axis: int,
            **kwargs: dict
        ) -> None:
        """
        Ungroup the GroupableMatrix along an axis by removing grouping metadata.

        Parameters
        ----------
        axis : int
            The axis along which values should be ungrouped.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def is_grouped(
            self, 
            axis: int, 
            **kwargs: dict
        ) -> bool:
        """
        Determine whether the Matrix has been sorted and grouped.

        Parameters
        ----------
        axis: int
            Axis along which to determine whether elements have been sorted and
            grouped.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        grouped : bool
            True or False indicating whether the Matrix has been sorted and
            grouped.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GroupableMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type GroupableMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GroupableMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,GroupableMatrix.__name__,type(v).__name__))
