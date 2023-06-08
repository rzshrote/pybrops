"""
Module defining interfaces and associated error checking routines for matrices
with taxa and trait metadata.
"""

__all__ = [
    "TaxaTraitMatrix",
    "check_is_TaxaTraitMatrix"
]

from pybrops.core.mat.TaxaMatrix import TaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix

class TaxaTraitMatrix(TaxaMatrix,TraitMatrix):
    """
    An abstract class for matrix wrapper objects with taxa and trait metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) TaxaMatrix
        2) TraitMatrix
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
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
def check_is_TaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type TaxaTraitMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, TaxaTraitMatrix):
        raise TypeError("variable '{0}' must be a of type '{1}' but received type '{2}'".format(vname,TaxaTraitMatrix.__name__,type(v).__name__))
