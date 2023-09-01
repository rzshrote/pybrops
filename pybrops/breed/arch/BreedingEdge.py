"""
Module defining interfaces and associated error checking routines for
breeding edges. Breeding edges define how information and germplasm flows
between breeding nodes.
"""

__all__ = [
    "BreedingEdge",
    "check_is_BreedingEdge"
]

from abc import ABCMeta


class BreedingEdge(metaclass=ABCMeta):
    """
    Abstract class defining interfaces for breeding edges. Breeding edges define
    how information and germplasm flows between breeding nodes.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################



################################## Utilities ###################################
def check_is_BreedingEdge(v: object, vname: str) -> None:
    """
    Check if object is of type BreedingEdge. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingEdge):
        raise TypeError("'%s' must be a BreedingEdge." % vname)
