"""
Module defining basal matrix interfaces and associated error checking routines
for phased haplotype matrices.
"""

__all__ = [
    "PhasedHaplotypeMatrix",
    "check_is_PhasedHaplotypeMatrix",
]

from abc import ABCMeta
from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.popgen.gmat.HaplotypeMatrix import HaplotypeMatrix

class PhasedHaplotypeMatrix(
        HaplotypeMatrix,
        PhasedTaxaVariantMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for phased genoypte matrix objects.

    The purpose of this abstract class is to merge the following interfaces:
        1) HaplotypeMatrix
        2) PhasedTaxaVariantMatrix
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_PhasedHaplotypeMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``PhasedHaplotypeMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, PhasedHaplotypeMatrix):
        raise TypeError("'%s' must be a PhasedHaplotypeMatrix." % vname)
