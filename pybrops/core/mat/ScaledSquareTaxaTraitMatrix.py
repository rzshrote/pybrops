"""
Module defining interfaces and associated error checking routines for scaled 
matrices with taxa axes that are square and a trait axis which is not square.
"""


from pybrops.core.mat.ScaledMatrix import ScaledMatrix
from pybrops.core.mat.SquareTaxaTraitMatrix import SquareTaxaTraitMatrix


class ScaledSquareTaxaTraitMatrix(
        SquareTaxaTraitMatrix,
        ScaledMatrix,
    ):
    """
    An abstract class for scaled matrix wrapper object with taxa axes which are
    square and a trait axis which is not square.

    The purpose of this abstract class is to merge the following interfaces:

        1. SquareTaxaTraitMatrix
        2. ScaledMatrix
    """
    pass



################################## Utilities ###################################
def check_is_ScaledSquareTaxaTraitMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``ScaledSquareTaxaTraitMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ScaledSquareTaxaTraitMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                ScaledSquareTaxaTraitMatrix.__name__,
                type(v).__name__
            )
        )
