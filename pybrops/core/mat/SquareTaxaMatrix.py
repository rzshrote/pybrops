"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square and with taxa metadata.
"""

__all__ = [
    "SquareTaxaMatrix",
    "check_is_SquareTaxaMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix

class SquareTaxaMatrix(
        SquareMatrix,
        TaxaMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper objects with square taxa metadata.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareMatrix
        2. TaxaMatrix
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## Square Metadata Properties ##############
    @property
    @abstractmethod
    def nsquare_taxa(self) -> int:
        """Number of taxa axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_taxa_axes(self) -> tuple:
        """Axis indices for taxa axes that are square."""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def square_taxa_axes_len(self) -> tuple:
        """Axis lengths for axes that are square."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    #################### Square Methods ####################
    @abstractmethod
    def is_square_taxa(
            self
        ) -> bool:
        """
        Determine whether the taxa axes lengths for the square axes are identical.

        Returns
        -------
        out : bool
            ``True`` if all square taxa axes are the same length.
            ``False`` if not all square taxa axes are the same length.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SquareTaxaMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTaxaMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTaxaMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareTaxaMatrix.__name__,type(v).__name__))
