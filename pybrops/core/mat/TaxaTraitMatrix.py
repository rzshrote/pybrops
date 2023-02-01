"""
Module defining interfaces and associated error checking routines for matrices
with taxa and trait metadata.
"""

from typing import Any
from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix

class TaxaTraitMatrix(TaxaMatrix,TraitMatrix):
    """
    An abstract class for matrix wrapper objects with taxa and trait metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) TaxaMatrix
        2) TraitMatrix
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict) -> None:
        """
        Constructor for the abstract class TaxaTraitMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for dependency injection.
        """
        super(TaxaTraitMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_TaxaTraitMatrix(v: Any) -> bool:
    """
    Determine whether an object is a TaxaTraitMatrix.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a TaxaTraitMatrix object instance.
    """
    return isinstance(v, TaxaTraitMatrix)

def check_is_TaxaTraitMatrix(v: Any, varname: str) -> None:
    """
    Check if object is of type TaxaTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TaxaTraitMatrix):
        raise TypeError("'{0}' must be a TaxaTraitMatrix".format(varname))
