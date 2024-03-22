"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square and with trait metadata.
"""

__all__ = [
    "SquareTraitMatrix",
    "check_is_SquareTraitMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix

class SquareTraitMatrix(
        SquareMatrix,
        TraitMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper objects with square axes and trait metadata.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareMatrix
        2. TraitMatrix
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @property
    @abstractmethod
    def nsquare_trait(self) -> int:
        """Number of trait axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_trait_axes(self) -> tuple:
        """Axis indices for trait axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_trait_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    @abstractmethod
    def is_square_trait(
            self
        ) -> bool:
        """
        Determine whether the trait axes lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square trait axes are the same length.
            ``False`` if not all square trait axes are the same length.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SquareTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareTraitMatrix.__name__,type(v).__name__))
