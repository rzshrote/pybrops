"""
Module for providing abstract interfaces and functions pertaining to square taxa
matrices.
"""

from pybrops.core.mat.SquareMatrix import SquareMatrix
from pybrops.core.mat.TaxaMatrix import TaxaMatrix

class SquareTaxaMatrix(SquareMatrix,TaxaMatrix):
    """Abstract class for SquareMatrix + TaxaMatrix fusion."""

    def __init__(self, **kwargs):
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
def is_SquareTaxaMatrix(obj):
    """
    Determine whether an object is a ``SquareTaxaMatrix``.

    Parameters
    ----------
    obj : object
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``SquareTaxaMatrix`` object instance.
    """
    return isinstance(obj, SquareTaxaMatrix)

def check_is_SquareTaxaMatrix(obj, objname):
    """
    Check if object is of type ``SquareTaxaMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(obj, SquareTaxaMatrix):
        raise TypeError("'{0}' must be a SquareTaxaMatrix".format(objname))

def cond_check_is_SquareTaxaMatrix(obj, objname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type ``SquareTaxaMatrix``. Otherwise raise
    ``TypeError``.

    Parameters
    ----------
    obj : object
        Any Python object to test.
    objname : str
        Name of variable to print in ``TypeError`` message.
    cond : function
        A function returning ``True`` or ``False`` for whether to test if ``obj``
        is a ``SquareTaxaMatrix``.
    """
    if cond(obj) and not isinstance(obj, SquareTaxaMatrix):
        raise TypeError("'{0}' must be a SquareTaxaMatrix".format(objname))
