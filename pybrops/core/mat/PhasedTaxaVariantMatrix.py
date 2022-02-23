"""
Module defining matrix interfaces and associated error checking routines for
matrices with phase, variant, and taxa metadata.
"""

from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybrops.core.mat.PhasedMatrix import PhasedMatrix

class PhasedTaxaVariantMatrix(TaxaVariantMatrix,PhasedMatrix):
    """
    An abstract class for matrix wrapper objects with phase, variant, and taxa
    metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) TaxaVariantMatrix
        2) PhasedMatrix
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class PhasedTaxaVariantMatrix.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
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
