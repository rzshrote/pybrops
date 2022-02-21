from pybrops.popgen.ldmat.LDMatrix import LDMatrix

# TODO: implementation
class R_LDMatrix(LDMatrix):
    """docstring for R_LDMatrix."""

    def __init__(self, arg):
        super(R_LDMatrix, self).__init__()
        self.arg = arg



################################################################################
################################## Utilities ###################################
################################################################################
def is_R_LDMatrix(v):
    return isinstance(v, R_LDMatrix)

def check_is_R_LDMatrix(v, varname):
    if not isinstance(v, R_LDMatrix):
        raise TypeError("'%s' must be a R_LDMatrix." % varname)

def cond_check_is_R_LDMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_R_LDMatrix(v, varname)
