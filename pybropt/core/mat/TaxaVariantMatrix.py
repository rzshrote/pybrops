from . import TaxaMatrix
from . import VariantMatrix

class TaxaVariantMatrix(TaxaMatrix,VariantMatrix):
    """Abstract class for TaxaMatrix + VariantMatrix class fusion."""

    def __init__(self, **kwargs):
        super(TaxaVariantMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_TaxaVariantMatrix(v):
    """
    Determine whether an object is a TaxaVariantMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TaxaVariantMatrix object instance.
    """
    return isinstance(v, TaxaVariantMatrix)

def check_is_TaxaVariantMatrix(v, varname):
    """
    Check if object is of type TaxaVariantMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_TaxaVariantMatrix(v):
        raise TypeError("'{0}' must be a TaxaVariantMatrix".format(varname))

def cond_check_is_TaxaVariantMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type TaxaVariantMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        TaxaVariantMatrix.
    """
    if cond(v):
        check_is_TaxaVariantMatrix(v, varname)