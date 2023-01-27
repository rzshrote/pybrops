"""
Module defining interfaces and associated error checking routines for
breeding edges. Breeding edges define how information and germplasm flows
between breeding nodes.
"""

from typing import Any


class BreedingEdge:
    """
    Abstract class defining interfaces for breeding edges. Breeding edges define
    how information and germplasm flows between breeding nodes.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        """
        Constructor for the abstract class BreedingEdge.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingEdge, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingEdge(v: Any) -> bool:
    """
    Determine whether an object is a BreedingEdge.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingEdge object instance.
    """
    return isinstance(v, BreedingEdge)

def check_is_BreedingEdge(v: Any, varname: str) -> None:
    """
    Check if object is of type BreedingEdge. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingEdge):
        raise TypeError("'%s' must be a BreedingEdge." % varname)
