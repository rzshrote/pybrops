"""
Module containing the abstract class GenotypingProtocol and its service functions.
"""

from abc import ABCMeta
from abc import abstractmethod
from typing import Optional

from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class GenotypingProtocol(
        metaclass = ABCMeta,
    ):
    """
    Abstract class defining genotyping protocols.

    The purpose of this abstract class is to define functionality for:
        1) Genotyping of individuals (converting a genome matrix to genotype matrix).
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def genotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> GenotypeMatrix:
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



################################## Utilities ###################################
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
