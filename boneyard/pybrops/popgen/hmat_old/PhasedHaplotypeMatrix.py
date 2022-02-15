from . import HaplotypeMatrix
from pybrops.core.mat import PhasedMatrix

class PhasedHaplotypeMatrix(HaplotypeMatrix,PhasedMatrix):
    """
    An abstract class for phased haplotype matrix objects.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PhasedHaplotypeMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhasedHaplotypeMatrix(v):
    return isinstance(v, PhasedHaplotypeMatrix)

def check_is_PhasedHaplotypeMatrix(v, varname):
    if not isinstance(v, PhasedHaplotypeMatrix):
        raise TypeError("'%s' must be a PhasedHaplotypeMatrix." % varname)

def cond_check_is_PhasedHaplotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PhasedHaplotypeMatrix(v, varname)
