"""
Module defining basal matrix interfaces and associated error checking routines
for phased genotype matrices.
"""

from pybrops.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

# NOTE: GenotypeMatrix MUST GO FIRST: causes method resolution error otherwise!
# NOTE: this is because
#       GenotypeMatrix extends TaxaVariantMatrix
#       PhasedTaxaVariantMatrix extends TaxaVariantMatrix
# I think a situation is caused where TaxaVariantMatrix is needed to be imported
# before itself due to its ranking in the MRO algorithm for multiple inheritance
class PhasedGenotypeMatrix(GenotypeMatrix,PhasedTaxaVariantMatrix):
    """
    An abstract class for phased genoypte matrix objects.

    The purpose of this abstract class is to merge the following interfaces:
        1) GenotypeMatrix
        2) PhasedTaxaVariantMatrix
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class PhasedGenotypeMatrix.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(PhasedGenotypeMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhasedGenotypeMatrix(v):
    """Return whether an object is a PhasedGenotypeMatrix or not"""
    return isinstance(v, PhasedGenotypeMatrix)

def check_is_PhasedGenotypeMatrix(v, varname):
    """Raise TypeError if object is not a PhasedGenotypeMatrix"""
    if not isinstance(v, PhasedGenotypeMatrix):
        raise TypeError("'%s' must be a PhasedGenotypeMatrix." % varname)

def cond_check_is_PhasedGenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a PhasedGenotypeMatrix"""
    if cond(v):
        check_is_PhasedGenotypeMatrix(v, varname)
