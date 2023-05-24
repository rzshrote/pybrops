"""
Module defining interfaces and associated error checking routines for matrices
with axes that are square and with taxa metadata.
"""

from typing import Any
from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix

class SquareTaxaMatrix(SquareMatrix,TaxaMatrix):
    """
    An abstract class for matrix wrapper objects with taxa and trait metadata.

    The purpose of this abstract class is to merge the following interfaces:
        1) SquareMatrix
        2) TaxaMatrix
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the SquareTaxaMatrix abstract class.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(SquareTaxaMatrix, self).__init__(**kwargs)



################################################################################
################################## Utilities ###################################
################################################################################
def is_SquareTaxaMatrix(v: object) -> bool:
    """
    Determine whether an object is a ``SquareTaxaMatrix``.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``SquareTaxaMatrix`` object instance.
    """
    return isinstance(v, SquareTaxaMatrix)

def check_is_SquareTaxaMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``SquareTaxaMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, SquareTaxaMatrix):
        raise TypeError("'{0}' must be a SquareTaxaMatrix".format(vname))
