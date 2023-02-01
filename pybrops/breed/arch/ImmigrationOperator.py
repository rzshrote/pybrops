"""
Module defining interfaces and associated error checking routines for
immigration operators.
"""

from typing import Any
from pybrops.breed.arch.BreedingEdge import BreedingEdge

class ImmigrationOperator(BreedingEdge):
    """
    Abstract class defining immigration operators.

    The purpose of this abstract class is to define functionality for:
        1) Protocols for immigration between different breeding nodes.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict) -> None:
        """
        Constructor for the abstract class ImmigrationOperator.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(ImmigrationOperator, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def immigrate(self, bnode, **kwargs: dict):
        """
        Immigrate individuals from a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object from which to pull individuals.

        Returns
        -------
        kwindiv : dict
            Dictionary of data from individuals selected.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_ImmigrationOperator(v: Any) -> bool:
    """
    Determine whether an object is a ImmigrationOperator.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a ImmigrationOperator object instance.
    """
    return isinstance(v, ImmigrationOperator)

def check_is_ImmigrationOperator(v: Any, varname: str) -> None:
    """
    Check if object is of type ImmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ImmigrationOperator):
        raise TypeError("'%s' must be a ImmigrationOperator." % varname)
