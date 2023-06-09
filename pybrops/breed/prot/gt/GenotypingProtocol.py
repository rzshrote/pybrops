"""
Module containing the abstract class GenotypingProtocol and its service functions.
"""

from typing import Any


class GenotypingProtocol:
    """
    Abstract class defining genotyping protocols.

    The purpose of this abstract class is to define functionality for:
        1) Genotyping of individuals (converting a genome matrix to genotype matrix).
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class GenotypingProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(GenotypingProtocol, self).__init__()

    ############################## Object Methods ##############################
    def genotype(self, pgmat, miscout, **kwargs: dict):
        """
        Genotype a genome.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A phased genotype matrix representing the whole simulated genome.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GenotypeMatrix
            A GenotypeMatrix of genotyped individuals.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_GenotypingProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type GenotypingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GenotypingProtocol):
        raise TypeError("'%s' must be a GenotypingProtocol." % vname)
