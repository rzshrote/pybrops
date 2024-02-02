"""
Module defining basal matrix interfaces and associated error checking routines
for phased genotype matrices.
"""

from abc import ABCMeta
from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

# NOTE: GenotypeMatrix MUST GO FIRST: causes method resolution error otherwise!
# NOTE: this is because
#       GenotypeMatrix extends TaxaVariantMatrix
#       PhasedTaxaVariantMatrix extends TaxaVariantMatrix
# I think a situation is caused where TaxaVariantMatrix is needed to be imported
# before itself due to its ranking in the MRO algorithm for multiple inheritance
class PhasedGenotypeMatrix(GenotypeMatrix,PhasedTaxaVariantMatrix,metaclass=ABCMeta):
    """
    An abstract class for phased genoypte matrix objects.

    The purpose of this abstract class is to merge the following interfaces:
        1) GenotypeMatrix
        2) PhasedTaxaVariantMatrix
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_PhasedGenotypeMatrix(v: object, vname: str) -> None:
    """Raise TypeError if object is not a PhasedGenotypeMatrix"""
    if not isinstance(v, PhasedGenotypeMatrix):
        raise TypeError("'%s' must be a PhasedGenotypeMatrix." % vname)
