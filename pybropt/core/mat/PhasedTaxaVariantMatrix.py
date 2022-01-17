from pybropt.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybropt.core.mat.PhasedMatrix import PhasedMatrix

class PhasedTaxaVariantMatrix(TaxaVariantMatrix,PhasedMatrix):
    """Abstract class for TaxaVariantMatrix + PhasedMatrix fusion."""

    def __init__(self, **kwargs):
        """
        Constructor for the abstract class PhasedTaxaVariantMatrix.
        """
        super(PhasedTaxaVariantMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhasedTaxaVariantMatrix(v):
    """
    Determine whether an object is a PhasedTaxaVariantMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PhasedTaxaVariantMatrix object instance.
    """
    return isinstance(v, PhasedTaxaVariantMatrix)

def check_is_PhasedTaxaVariantMatrix(v, varname):
    """
    Check if object is of type PhasedTaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_PhasedTaxaVariantMatrix(v):
        raise TypeError("'{0}' must be a PhasedTaxaVariantMatrix".format(varname))

def cond_check_is_PhasedTaxaVariantMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type PhasedTaxaVariantMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        PhasedTaxaVariantMatrix.
    """
    if cond(v):
        check_is_PhasedTaxaVariantMatrix(v, varname)
