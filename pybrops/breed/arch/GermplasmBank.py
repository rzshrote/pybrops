"""
Module defining interfaces and associated error checking routines for
germplasm bank representation.
"""

from typing import Any
from pybrops.breed.arch.BreedingNode import BreedingNode

class GermplasmBank(BreedingNode):
    """
    Abstract class defining a germplasm bank.

    The purpose of this abstract class is to define functionality for:
        1) Container storage for the germplasm bank.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the abstract class GermplasmBank.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(GermplasmBank, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_GermplasmBank(v: Any) -> bool:
    """
    Determine whether an object is a GermplasmBank.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GermplasmBank object instance.
    """
    return isinstance(v, GermplasmBank)

def check_is_GermplasmBank(v: Any, varname: str) -> None:
    """
    Check if object is of type GermplasmBank. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GermplasmBank):
        raise TypeError("'%s' must be a GermplasmBank." % varname)
