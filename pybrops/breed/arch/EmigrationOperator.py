"""
Module defining interfaces and associated error checking routines for
emigration operators.
"""

__all__ = [
    "EmigrationOperator",
    "check_is_EmigrationOperator"
]

from abc import ABCMeta, abstractmethod
from pybrops.breed.arch.BreedingEdge import BreedingEdge
from pybrops.breed.arch.BreedingNode import BreedingNode

class EmigrationOperator(BreedingEdge,metaclass=ABCMeta):
    """
    Abstract class defining immigration operators.

    The purpose of this abstract class is to define functionality for:
        1) Protocols for emigration between different breeding nodes.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    @abstractmethod
    def emigrate(
            self, 
            bnode: BreedingNode, 
            **kwargs: dict
        ) -> None:
        """
        Emigrate individuals to a BreedingNode.

        Parameters
        ----------
        bnode : BreedingNode
            A BreedingNode object to add individuals.
        kwargs : dict
            Dictionary of data for individuals to be added.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_EmigrationOperator(v: object, vname: str) -> None:
    """
    Check if object is of type EmigrationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, EmigrationOperator):
        raise TypeError("'%s' must be a EmigrationOperator." % vname)
