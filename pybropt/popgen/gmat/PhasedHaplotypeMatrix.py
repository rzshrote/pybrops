from pybropt.core.mat.PhasedTaxaVariantMatrix import PhasedTaxaVariantMatrix
from pybropt.popgen.gmat.HaplotypeMatrix import HaplotypeMatrix

class PhasedHaplotypeMatrix(PhasedTaxaVariantMatrix,HaplotypeMatrix):
    """docstring for PhasedHaplotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class PhasedHaplotypeMatrix.
        """
        super(PhasedHaplotypeMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhasedHaplotypeMatrix(v):
    """Return whether an object is a PhasedHaplotypeMatrix or not"""
    return isinstance(v, PhasedHaplotypeMatrix)

def check_is_PhasedHaplotypeMatrix(v, varname):
    """Raise TypeError if object is not a PhasedHaplotypeMatrix"""
    if not isinstance(v, PhasedHaplotypeMatrix):
        raise TypeError("'%s' must be a PhasedHaplotypeMatrix." % varname)

def cond_check_is_PhasedHaplotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a PhasedHaplotypeMatrix"""
    if cond(v):
        check_is_PhasedHaplotypeMatrix(v, varname)
