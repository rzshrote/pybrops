from . import GenotypeMatrix
from pybropt.core.mat import TaxaMatrix
from pybropt.core.mat import VariantMatrix
from pybropt.core.mat import PrunableMatrix

class GenotypeVariantMatrix(GenotypeMatrix,TaxaMatrix,VariantMatrix,PrunableMatrix):
    """docstring for GenotypeVariantMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        GenotypeVariantMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GenotypeVariantMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenotypeVariantMatrix(v):
    return isinstance(v, GenotypeVariantMatrix)

def check_is_GenotypeVariantMatrix(v, varname):
    if not isinstance(v, GenotypeVariantMatrix):
        raise TypeError("'%s' must be a GenotypeVariantMatrix." % varname)

def cond_check_is_GenotypeVariantMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenotypeVariantMatrix(v, varname)
