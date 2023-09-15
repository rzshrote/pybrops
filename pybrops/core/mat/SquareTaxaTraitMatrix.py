"""
Module defining interfaces and associated error checking routines for matrices
with taxa axes that are square and a trait axis which is not square.
"""

__all__ = [
    "SquareTaxaTraitMatrix",
    "check_is_SquareTaxaTraitMatrix",
]

from abc import ABCMeta
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix


class SquareTaxaTraitMatrix(SquareTaxaMatrix,TraitMatrix,metaclass=ABCMeta):
    """
    An abstract class for matrix wrapper object with taxa axes which are square
    and a trait axis which is not square.

    The purpose of this abstract class is to merge the following interfaces:
        1) SquareTaxaMatrix
        2) TraitMatrix
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_SquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTaxaTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTaxaTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareTaxaTraitMatrix.__name__,type(v).__name__))
