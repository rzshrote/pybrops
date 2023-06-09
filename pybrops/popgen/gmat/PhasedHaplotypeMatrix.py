"""
Module defining basal matrix interfaces and associated error checking routines
for phased haplotype matrices.
"""

__all__ = [
    "PhasedHaplotypeMatrix",
    "check_is_PhasedHaplotypeMatrix"
]

from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.popgen.gmat.HaplotypeMatrix import HaplotypeMatrix

class PhasedHaplotypeMatrix(HaplotypeMatrix,PhasedTaxaVariantMatrix):
    """
    An abstract class for phased genoypte matrix objects.

    The purpose of this abstract class is to merge the following interfaces:
        1) HaplotypeMatrix
        2) PhasedTaxaVariantMatrix
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class PhasedHaplotypeMatrix.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(PhasedHaplotypeMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_PhasedHaplotypeMatrix(v: object, vname: str) -> None:
    """Raise TypeError if object is not a PhasedHaplotypeMatrix"""
    if not isinstance(v, PhasedHaplotypeMatrix):
        raise TypeError("'%s' must be a PhasedHaplotypeMatrix." % vname)
