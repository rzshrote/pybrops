from . import MutableMatrix
from . import SortableMatrixInterface

class SortableMatrix(MutableMatrix,SortableMatrixInterface):
    """docstring for SortableMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        SortableMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(SortableMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_SortableMatrix(v):
    return isinstance(v, SortableMatrix)

def check_is_SortableMatrix(v, varname):
    if not isinstance(v, SortableMatrix):
        raise TypeError("'%s' must be a SortableMatrix." % varname)

def cond_check_is_SortableMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SortableMatrix(v, varname)
