"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square and with taxa metadata.
"""

__all__ = [
    "SquareTaxaMatrix",
    "check_is_SquareTaxaMatrix",
]

from abc import ABCMeta
from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix

class SquareTaxaMatrix(SquareMatrix,TaxaMatrix,metaclass=ABCMeta):
    """
    An abstract class for matrix wrapper objects with taxa and trait metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) SquareMatrix
        2) TaxaMatrix
    """

    ########################## Special Object Methods ##########################



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
