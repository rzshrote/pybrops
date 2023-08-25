"""
Module defining matrix interfaces and associated error checking routines for
matrices with phase, variant, and taxa metadata.
"""

__all__ = [
    "PhasedTaxaVariantMatrix",
    "check_is_PhasedTaxaVariantMatrix"
]

from abc import ABCMeta
from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybrops.core.mat.PhasedMatrix import PhasedMatrix

class PhasedTaxaVariantMatrix(TaxaVariantMatrix,PhasedMatrix,metaclass=ABCMeta):
    """
    An abstract class for matrix wrapper objects with phase, variant, and taxa
    metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) TaxaVariantMatrix
        2) PhasedMatrix
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_PhasedTaxaVariantMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type PhasedTaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PhasedTaxaVariantMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,PhasedTaxaVariantMatrix.__name__,type(v).__name__))
