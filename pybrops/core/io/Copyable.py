"""
Module defining an interface for objects which are copyable and associated error checking routines.
"""

__all__ = [
    "Copyable",
    "check_is_Copyable",
]

from abc import ABCMeta, abstractmethod
from typing import Optional


class Copyable(
        metaclass = ABCMeta,
    ):
    """
    Abstract class for defining objects which can be copied.

    This abstract class defines four functions with the following purposes:

    - ``__copy__`` and ``__deepcopy__`` implement the Python ``copy`` interface.
    - ``copy`` and ``deepcopy`` allow for the direct copying of objects.
    """

    ########################## Special Object Methods ##########################
    @abstractmethod
    def __copy__(
            self
        ) -> 'Copyable':
        """
        Make a shallow copy of the ``Copyable`` object.

        Returns
        -------
        out : Copyable
            A shallow copy of the ``Copyable`` object.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __deepcopy__(
            self,
            memo: Optional[dict],
        ) -> 'Copyable':
        """
        Make a deep copy of the ``Copyable`` object.

        Parameters
        ----------
        memo : dict, None
            An optional dictionary of memo metadata.

        Returns
        -------
        out : Copyable
            A deep copy of the ``Copyable`` object.
        """
        raise NotImplementedError("method is abstract")

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    @abstractmethod
    def copy(
            self
        ) -> 'Copyable':
        """
        Make a shallow copy of the ``Copyable`` object.

        Returns
        -------
        out : Copyable
            A shallow copy of the ``Copyable`` object.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def deepcopy(
            self,
            memo: Optional[dict]
        ) -> 'Copyable':
        """
        Make a deep copy of the ``Copyable`` object.

        Parameters
        ----------
        memo : dict, None
            An optional dictionary of memo metadata.

        Returns
        -------
        out : Copyable
            A deep copy of the ``Copyable`` object.
        """
        raise NotImplementedError("method is abstract")

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_Copyable(v: object, vname: str) -> None:
    """
    Check if object is of type ``Copyable``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, Copyable):
        raise TypeError(
            "variable ``{0}`` must be of type ``Copyable`` but received type ``{1}``".format(
                vname,
                type(v).__name__,
            )
        )
