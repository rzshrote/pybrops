"""
Module defining interfaces and associated error checking routines for matrices
with taxa and variant metadata.
"""

__all__ = [
    "TaxaVariantMatrix",
    "check_is_TaxaVariantMatrix",
]

from abc import ABCMeta
from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.VariantMatrix import VariantMatrix

class TaxaVariantMatrix(
        TaxaMatrix,
        VariantMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrix wrapper objects with taxa and variant metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) TaxaMatrix
        2) VariantMatrix
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_TaxaVariantMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type TaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TaxaVariantMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,TaxaVariantMatrix.__name__,type(v).__name__))
