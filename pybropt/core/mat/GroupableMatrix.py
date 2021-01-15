from . import SortableMatrix
from . import GroupableMatrixInterface

class GroupableMatrix(SortableMatrix,GroupableMatrixInterface):
    """docstring for GroupableMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GroupableMatrix, self).__init__()




################################################################################
################################## Utilities ###################################
################################################################################
def is_GroupableMatrix(v):
    return isinstance(v, GroupableMatrix)

def check_is_GroupableMatrix(v, varname):
    if not isinstance(v, GroupableMatrix):
        raise TypeError("'%s' must be a GroupableMatrix." % varname)

def cond_check_is_GroupableMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GroupableMatrix(v, varname)
