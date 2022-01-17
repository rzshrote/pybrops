from pybropt.popgen.ldmat.LDMatrix import LDMatrix

# TODO: implementation
class Dprime_LDMatrix(LDMatrix):
    """docstring for Dprime_LDMatrix."""

    def __init__(self, arg):
        super(Dprime_LDMatrix, self).__init__()
        self.arg = arg



################################################################################
################################## Utilities ###################################
################################################################################
def is_Dprime_LDMatrix(v):
    return isinstance(v, Dprime_LDMatrix)

def check_is_Dprime_LDMatrix(v, varname):
    if not isinstance(v, Dprime_LDMatrix):
        raise TypeError("'%s' must be a Dprime_LDMatrix." % varname)

def cond_check_is_Dprime_LDMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_Dprime_LDMatrix(v, varname)
