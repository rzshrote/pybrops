"""
Module defining interfaces and associated error checking routines for breeding
graphs. Breeding program graphs represent breeding programs with multiple sub-
populations. Germplasm and information is passed between breeding nodes through
breeding edges.
"""

__all__ = [
    "BreedingGraph",
    "check_is_BreedingGraph"
]

from typing import Any

class BreedingGraph:
    """
    Abstract class defining interfaces for breeding graphs. Breeding graphs
    represent breeding programs with multiple subpopulations. Germplasm and
    information is passed between breeding nodes through breeding edges.

    The purpose of this abstract class is to provide functionality for:
        1) Graph representation of the entire breeding graph.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class BreedingGraph.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingGraph, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Graph properties ###################
    @property
    def graph(self) -> Any:
        """Graph data structure."""
        raise NotImplementedError("property is abstract")
    @graph.setter
    def graph(self, value: Any) -> None:
        """Set graph data structure."""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################## Utilities ###################################
def check_is_BreedingGraph(v: object, vname: str) -> None:
    """
    Check if object is of type BreedingGraph. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingGraph):
        raise TypeError("'%s' must be a BreedingGraph." % vname)
