"""
Module defining interfaces and associated error checking routines for matrices 
with square taxa axes and square trait axes.
"""

__all__ = [
    "SquareTaxaSquareTraitMatrix",
    "check_is_SquareTaxaSquareTraitMatrix",
]

from abc import ABCMeta

from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.SquareTraitMatrix import SquareTraitMatrix

class SquareTaxaSquareTraitMatrix(
        SquareTaxaMatrix,
        SquareTraitMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper object with square taxa and square 
    trait metadata.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareTaxaMatrix
        2. SquareTraitMatrix
    """
    pass



################################## Utilities ###################################
def check_is_SquareTaxaSquareTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTaxaSquareTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTaxaSquareTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,SquareTaxaSquareTraitMatrix.__name__,type(v).__name__))
