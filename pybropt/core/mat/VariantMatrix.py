from . import GroupableMatrix
from . import VariantMatrixInterface

class VariantMatrix(GroupableMatrix,VariantMatrixInterface):
    """docstring for VariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        VariantMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(VariantMatrix, self).__init__(**kwargs)




################################################################################
################################## Utilities ###################################
################################################################################
def is_VariantMatrix(v):
    return isinstance(v, VariantMatrix)

def check_is_VariantMatrix(v, varname):
    if not isinstance(v, VariantMatrix):
        raise TypeError("'%s' must be a VariantMatrix." % varname)

def cond_check_is_VariantMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_VariantMatrix(v, varname)
