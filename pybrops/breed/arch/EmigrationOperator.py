"""
Module defining interfaces and associated error checking routines for
emigration operators.
"""

__all__ = [
    "EmigrationOperator",
    "check_is_EmigrationOperator"
]

from pybrops.breed.arch.BreedingEdge import BreedingEdge

class EmigrationOperator(BreedingEdge):
    """
    Abstract class defining immigration operators.

    The purpose of this abstract class is to define functionality for:
        1) Protocols for emigration between different breeding nodes.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class EmigrationOperator.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(EmigrationOperator, self).__init__(**kwargs)

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    def emigrate(self, bnode, **kwargs: dict):
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
