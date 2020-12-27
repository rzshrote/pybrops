from . import LDMatrix

class D_LDMatrix(LDMatrix):
    """docstring for D_LDMatrix."""

    def __init__(self, arg):
        super(D_LDMatrix, self).__init__()
        self.arg = arg



################################################################################
################################## Utilities ###################################
################################################################################
def is_D_LDMatrix(v):
    return isinstance(v, D_LDMatrix)

def check_is_D_LDMatrix(v, varname):
    if not isinstance(v, D_LDMatrix):
        raise TypeError("'%s' must be a D_LDMatrix." % varname)

def cond_check_is_D_LDMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_D_LDMatrix(v, varname)
