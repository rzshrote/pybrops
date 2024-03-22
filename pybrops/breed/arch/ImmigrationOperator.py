"""
Module defining interfaces and associated error checking routines for
immigration operators.
"""

__all__ = [
    "ImmigrationOperator",
    "check_is_ImmigrationOperator",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.breed.arch.BreedingEdge import BreedingEdge
from pybrops.breed.arch.BreedingNode import BreedingNode

class ImmigrationOperator(BreedingEdge,metaclass=ABCMeta):
    """
    Abstract class defining immigration operators.

    The purpose of this abstract class is to define functionality for:
        1) Protocols for immigration between different breeding nodes.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    @abstractmethod
    def immigrate(
            self, 
            bnode: BreedingNode, 
            **kwargs: dict
        ) -> dict:
        """
        Immigrate individuals from a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object from which to pull individuals.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        kwindiv : dict
            Dictionary of data from individuals selected.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_ImmigrationOperator(v: object, vname: str) -> None:
    """
    Check if object is of type ImmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, ImmigrationOperator):
        raise TypeError("'%s' must be a ImmigrationOperator." % vname)
