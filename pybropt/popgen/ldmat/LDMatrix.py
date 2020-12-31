from pybropt.core.mat import Matrix

class LDMatrix(Matrix):
    """docstring for LDMatrix."""

    def __init__(self, **kwargs):
        super(LDMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_LDMatrix(v):
    return isinstance(v, LDMatrix)

def check_is_LDMatrix(v, varname):
    if not isinstance(v, LDMatrix):
        raise TypeError("'%s' must be a LDMatrix." % varname)

def cond_check_is_LDMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_LDMatrix(v, varname)
