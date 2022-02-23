"""
Module defining interfaces and associated error checking routines for matrices
with taxa and variant metadata.
"""

from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.VariantMatrix import VariantMatrix

class TaxaVariantMatrix(TaxaMatrix,VariantMatrix):
    """
    An abstract class for matrix wrapper objects with taxa and variant metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) TaxaMatrix
        2) VariantMatrix
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class TaxaVariantMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for dependency injection.
        """
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
