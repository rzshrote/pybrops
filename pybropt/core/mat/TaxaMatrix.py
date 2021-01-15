from . import GroupableMatrix
from . import TaxaMatrixInterface

class TaxaMatrix(GroupableMatrix,TaxaMatrixInterface):
    """docstring for TaxaMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        TaxaMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(TaxaMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_TaxaMatrix(v):
    return isinstance(v, TaxaMatrix)

def check_is_TaxaMatrix(v, varname):
    if not isinstance(v, TaxaMatrix):
        raise TypeError("'%s' must be a TaxaMatrix." % varname)

def cond_check_is_TaxaMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_TaxaMatrix(v, varname)
