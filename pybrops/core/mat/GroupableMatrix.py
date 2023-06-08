"""
Module defining interfaces and associated error checking routines for matrices
that can have their axes grouped.
"""

__all__ = [
    "GroupableMatrix",
    "check_is_GroupableMatrix"
]

from pybrops.core.mat.SortableMatrix import SortableMatrix

class GroupableMatrix(SortableMatrix):
    """
    An abstract class for groupable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix axis grouping routines.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        GroupableMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GroupableMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Grouping Methods ###################
    def group(
            self, 
            axis: int, 
            **kwargs: dict
        ) -> None:
        """
        Sort the Matrix, then populate grouping indices.

        Parameters
        ----------
        axis : int
            The axis along which values are grouped.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

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



################################################################################
################################## Utilities ###################################
################################################################################
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
        raise TypeError("'{0}' must be a GroupableMatrix".format(vname))
