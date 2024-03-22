"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square.
"""

__all__ = [
    "SquareMatrix",
    "check_is_SquareMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.Matrix import Matrix

class SquareMatrix(Matrix,metaclass=ABCMeta):
    """
    An abstract class for square matrices. A "square matrix" is defined as a
    matrix that has the same axis metadata associated with two or more axes.
    For example::

        This is a square matrix since metadata applies to axes 0 and 1:
               taxa
             +-------+
        taxa | (n,n) |
             +-------+

        This is not a square matrix since unique metadata applies to each axis:
               vrnt
             +-------+
        taxa | (n,p) |  Where: n == p
             +-------+

    The purpose of this abstract class is to provide base functionality for:

        1. Square matrix axis metadata.
        2. Determination of square matrix conformity.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @property
    @abstractmethod
    def nsquare(self) -> int:
        """Number of axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_axes(self) -> tuple:
        """Axis indices for axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    @abstractmethod
    def is_square(
            self
        ) -> bool:
        """
        Determine whether the axes lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square axes are the same length.
            ``False`` if not all square axes are the same length.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SquareMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareMatrix.__name__,type(v).__name__))
