from . import Matrix
from . import MutableMatrixInterface

class MutableMatrix(Matrix,MutableMatrixInterface):
    """docstring for MutableMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(MutableMatrix, self).__init__(**kwargs)


################################################################################
################################## Utilities ###################################
################################################################################
def is_MutableMatrix(v):
    return isinstance(v, MutableMatrix)

def check_is_MutableMatrix(v, varname):
    if not isinstance(v, MutableMatrix):
        raise TypeError("'%s' must be a MutableMatrix." % varname)

def cond_check_is_MutableMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_MutableMatrix(v, varname)
