"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square.
"""

from typing import Any
from pybrops.core.mat.Matrix import Matrix

class SquareMatrix(Matrix):
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
        1) Square matrix axis metadata.
        2) Determination of square matrix conformity.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the SquareMatrix abstract class.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(SquareMatrix, self).__init__(**kwargs)

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @property
    def nsquare(self) -> int:
        """Number of axes that are square."""
        raise NotImplementedError("property is abstract")
    @nsquare.setter
    def nsquare(self, value: int) -> None:
        """Set the number of axes that are square"""
        raise NotImplementedError("property is abstract")
    @nsquare.deleter
    def nsquare(self) -> None:
        """Delete the number of axes that are square"""
        raise NotImplementedError("property is abstract")
    
    @property
    def square_axes(self) -> tuple:
        """Axis indices for axes that are square."""
        raise NotImplementedError("property is abstract")
    @square_axes.setter
    def square_axes(self, value: tuple) -> None:
        """Set axis indices for axes that are square"""
        raise NotImplementedError("property is abstract")
    @square_axes.deleter
    def square_axes(self) -> None:
        """Delete axis indices for axes that are square"""
        raise NotImplementedError("property is abstract")
    
    @property
    def square_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        raise NotImplementedError("property is abstract")
    @square_axes_len.setter
    def square_axes_len(self, value: tuple) -> None:
        """Set axis lengths for axes that are square"""
        raise NotImplementedError("property is abstract")
    @square_axes_len.deleter
    def square_axes_len(self) -> None:
        """Delete axis lengths for axes that are square"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    def is_square(
            self
        ) -> bool:
        """
        Determine whether the axis lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square axes are the same length.
            ``False`` if not all square axes are the same length.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_SquareMatrix(v: object) -> bool:
    """
    Determine whether an object is a ``SquareMatrix``.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``SquareMatrix`` object instance.
    """
    return isinstance(v, SquareMatrix)

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
        raise TypeError("'{0}' must be a SquareMatrix".format(vname))
