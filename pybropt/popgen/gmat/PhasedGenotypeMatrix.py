from pybropt.core.mat import PhasedTaxaVariantMatrix
from . import GenotypeMatrix

# NOTE: GenotypeMatrix MUST GO FIRST: causes method resolution error otherwise!
# NOTE: this is because
#       GenotypeMatrix extends TaxaVariantMatrix
#       PhasedTaxaVariantMatrix extends TaxaVariantMatrix
# I think a situation is caused where TaxaVariantMatrix is needed to be imported
# before itself due to its ranking in the MRO algorithm for multiple inheritance
class PhasedGenotypeMatrix(GenotypeMatrix,PhasedTaxaVariantMatrix):
    """docstring for PhasedGenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PhasedGenotypeMatrix, self).__init__(**kwargs)

################################################################################
################################## Utilities ###################################
################################################################################
def is_PhasedGenotypeMatrix(v):
    return isinstance(v, PhasedGenotypeMatrix)

def check_is_PhasedGenotypeMatrix(v, varname):
    if not isinstance(v, PhasedGenotypeMatrix):
        raise TypeError("'%s' must be a PhasedGenotypeMatrix." % varname)

def cond_check_is_PhasedGenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PhasedGenotypeMatrix(v, varname)
