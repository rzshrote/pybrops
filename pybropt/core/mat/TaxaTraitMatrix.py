from . import TaxaMatrix
from . import TraitMatrix

class TaxaTraitMatrix(TaxaMatrix,TraitMatrix):
    """Abstract class for TaxaMatrix + TraitMatrix fusion."""

    def __init__(self, **kwargs):
        super(TaxaTraitMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_TaxaTraitMatrix(v):
    """
    Determine whether an object is a TaxaTraitMatrix.

    Parameters
    ----------
    v : any object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TaxaTraitMatrix object instance.
    """
    return isinstance(v, TaxaTraitMatrix)

def check_is_TaxaTraitMatrix(v, varname):
    """
    Check if object is of type TaxaTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not is_TaxaTraitMatrix(v):
        raise TypeError("'{0}' must be a TaxaTraitMatrix".format(varname))

def cond_check_is_TaxaTraitMatrix(v, varname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type TaxaTraitMatrix. Otherwise raise
    TypeError.

    Parameters
    ----------
    v : any object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a
        TaxaTraitMatrix.
    """
    if cond(v):
        check_is_TaxaTraitMatrix(v, varname)
