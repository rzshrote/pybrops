"""
Module defining interfaces and associated error checking routines for
germplasm bank representation.
"""

__all__ = [
    "GermplasmBank",
    "check_is_GermplasmBank"
]

from abc import ABCMeta
from pybrops.breed.arch.BreedingNode import BreedingNode

class GermplasmBank(BreedingNode,metaclass=ABCMeta):
    """
    Abstract class defining a germplasm bank.

    The purpose of this abstract class is to define functionality for:
        1) Container storage for the germplasm bank.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################



################################## Utilities ###################################
def check_is_GermplasmBank(v: object, vname: str) -> None:
    """
    Check if object is of type GermplasmBank. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GermplasmBank):
        raise TypeError("'%s' must be a GermplasmBank." % vname)
