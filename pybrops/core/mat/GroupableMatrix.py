"""
Module defining interfaces and associated error checking routines for matrices
that can have their axes grouped.
"""

__all__ = [
    "GroupableMatrix",
    "check_is_GroupableMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.SortableMatrix import SortableMatrix

class GroupableMatrix(
        SortableMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for groupable matrix wrapper objects.

    The purpose of this abstract class is to provide base functionality for:
        1) Matrix in-place matrix axis grouping routines.
    """

    ########################## Special Object Methods ##########################

    ############## Forward numeric operators ###############
    ### __add__         inherited from ``SortableMatrix``
    ### __sub__         inherited from ``SortableMatrix``
    ### __mul__         inherited from ``SortableMatrix``
    ### __matmul__      inherited from ``SortableMatrix``
    ### __truediv__     inherited from ``SortableMatrix``
    ### __floordiv__    inherited from ``SortableMatrix``
    ### __mod__         inherited from ``SortableMatrix``
    ### __divmod__      inherited from ``SortableMatrix``
    ### __pow__         inherited from ``SortableMatrix``
    ### __lshift__      inherited from ``SortableMatrix``
    ### __rshift__      inherited from ``SortableMatrix``
    ### __and__         inherited from ``SortableMatrix``
    ### __xor__         inherited from ``SortableMatrix``
    ### __or__          inherited from ``SortableMatrix``

    ############# Backwards numeric operators ##############
    ### __radd__        inherited from ``SortableMatrix``
    ### __rsub__        inherited from ``SortableMatrix``
    ### __rmul__        inherited from ``SortableMatrix``
    ### __rmatmul__     inherited from ``SortableMatrix``
    ### __rtruediv__    inherited from ``SortableMatrix``
    ### __rfloordiv__   inherited from ``SortableMatrix``
    ### __rmod__        inherited from ``SortableMatrix``
    ### __rdivmod__     inherited from ``SortableMatrix``
    ### __rlshift__     inherited from ``SortableMatrix``
    ### __rrshift__     inherited from ``SortableMatrix``
    ### __rand__        inherited from ``SortableMatrix``
    ### __rxor__        inherited from ``SortableMatrix``
    ### __ror__         inherited from ``SortableMatrix``

    ############# Augmented numeric operators ##############
    ### __iadd__        inherited from ``SortableMatrix``
    ### __isub__        inherited from ``SortableMatrix``
    ### __imul__        inherited from ``SortableMatrix``
    ### __imatmul__     inherited from ``SortableMatrix``
    ### __itruediv__    inherited from ``SortableMatrix``
    ### __ifloordiv__   inherited from ``SortableMatrix``
    ### __imod__        inherited from ``SortableMatrix``
    ### __ipow__        inherited from ``SortableMatrix``
    ### __ilshift__     inherited from ``SortableMatrix``
    ### __irshift__     inherited from ``SortableMatrix``
    ### __iand__        inherited from ``SortableMatrix``
    ### __ixor__        inherited from ``SortableMatrix``
    ### __ior__         inherited from ``SortableMatrix``

    ################## Logical operators ###################
    ### __lt__          inherited from ``SortableMatrix``
    ### __le__          inherited from ``SortableMatrix``
    ### __eq__          inherited from ``SortableMatrix``
    ### __ne__          inherited from ``SortableMatrix``
    ### __gt__          inherited from ``SortableMatrix``
    ### __ge__          inherited from ``SortableMatrix``

    ################# Container operators ##################
    ### __len__         inherited from ``SortableMatrix``
    ### __getitem__     inherited from ``SortableMatrix``
    ### __setitem__     inherited from ``SortableMatrix``
    ### __delitem__     inherited from ``SortableMatrix``
    ### __iter__        inherited from ``SortableMatrix``

    #################### Matrix copying ####################
    ### __copy__        inherited from ``SortableMatrix``
    ### __deepcopy__    inherited from ``SortableMatrix``

    ########### Miscellaneous special functions ############
    ### __repr__        inherited from ``SortableMatrix``

    ############################ Object Properties #############################

    ################## Matrix Properties ###################
    ### mat             inherited from ``SortableMatrix``

    ############## Matrix Metadata Properties ##############
    ### mat_ndim        inherited from ``SortableMatrix``
    ### mat_shape       inherited from ``SortableMatrix``

    ############################## Object Methods ##############################

    #################### Matrix copying ####################
    ### copy            inherited from ``SortableMatrix``
    ### deepcopy        inherited from ``SortableMatrix``

    ######### Matrix element copy-on-manipulation ##########
    ### adjoin          inherited from ``SortableMatrix``
    ### delete          inherited from ``SortableMatrix``
    ### insert          inherited from ``SortableMatrix``
    ### select          inherited from ``SortableMatrix``

    ######### Matrix element in-place-manipulation #########
    ### append          inherited from ``SortableMatrix``
    ### remove          inherited from ``SortableMatrix``
    ### incorp          inherited from ``SortableMatrix``

    ################### Sorting Methods ####################
    ### lexsort         inherited from ``SortableMatrix``
    ### reorder         inherited from ``SortableMatrix``
    ### sort            inherited from ``SortableMatrix``

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

    ############################## Class Methods ###############################

    ######### Matrix element copy-on-manipulation ##########
    ### concat          inherited from ``SortableMatrix``



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
