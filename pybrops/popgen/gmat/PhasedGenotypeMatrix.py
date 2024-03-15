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
class PhasedGenotypeMatrix(
        GenotypeMatrix,
        PhasedTaxaVariantMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for phased genoypte matrix objects.

    The purpose of this abstract class is to merge the following interfaces:
        1) GenotypeMatrix
        2) PhasedTaxaVariantMatrix
    """

    ########################## Special Object Methods ##########################



################################## Utilities ###################################
def check_is_PhasedGenotypeMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``PhasedGenotypeMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, PhasedGenotypeMatrix):
        raise TypeError("'%s' must be a PhasedGenotypeMatrix." % vname)

def check_PhasedGenotypeMatrix_has_vrnt_mask(v: PhasedGenotypeMatrix, vname: str) -> None:
    """
    Check if a PhasedGenotypeMatrix has the vrnt_mask field that is non-None.

    Parameters
    ----------
    v : PhasedGenotypeMatrix
        A phased genotype matrix.
    vname : str
        Name of the variable to print in ``ValueError`` message.
    """
    if v.vrnt_mask is None:
        raise ValueError(
            "{0} '{0}' must have 'vrnt_mask' (variant mask) assigned".format(
                type(v).__name__,
                vname
            )
        )