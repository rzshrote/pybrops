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
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
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
def is_BreedingEdge(v: object) -> bool:
    """
    Determine whether an object is a BreedingEdge.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingEdge object instance.
    """
    return isinstance(v, BreedingEdge)

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
        raise TypeError("'%s' must be a BreedingEdge." % varname)
