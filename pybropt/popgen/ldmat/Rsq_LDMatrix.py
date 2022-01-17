from pybropt.popgen.ldmat.LDMatrix import LDMatrix

# TODO: implementation
class Rsq_LDMatrix(LDMatrix):
    """docstring for Rsq_LDMatrix."""

    def __init__(self, arg):
        super(Rsq_LDMatrix, self).__init__()
        self.arg = arg



################################################################################
################################## Utilities ###################################
################################################################################
def is_Rsq_LDMatrix(v):
    return isinstance(v, Rsq_LDMatrix)

def check_is_Rsq_LDMatrix(v, varname):
    if not isinstance(v, Rsq_LDMatrix):
        raise TypeError("'%s' must be a Rsq_LDMatrix." % varname)

def cond_check_is_Rsq_LDMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_Rsq_LDMatrix(v, varname)
