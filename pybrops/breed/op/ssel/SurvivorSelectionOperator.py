"""
Module defining interfaces and associated error checking routines for
breeding program survivor selection operators.
"""

__all__ = [
    "SurvivorSelectionOperator",
    "check_is_SurvivorSelectionOperator",
]

from abc import ABCMeta, abstractmethod
from typing import Optional


class SurvivorSelectionOperator(
        metaclass = ABCMeta,
    ):
    """
    Abstract class defining interfaces for survivor selection within an entire
    breeding program.

    The purpose of this abstract class is to provide functionality for:
        1) Survivor selection for an entire breeding program.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def sselect(
            self, 
            genome: dict, 
            geno: dict, 
            pheno: dict, 
            bval: dict, 
            gmod: dict, 
            t_cur: int, 
            t_max: int, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> tuple:
        """
        Select progeny survivors in a breeding program.

        Parameters
        ----------
        genome : dict
            Dictionary of genomes for the breeding program.
        geno : dict
            Dictionary of genotypes for the breeding program.
        pheno : dict
            Dictionary of phenotypes for the breeding program.
        bval : dict
            Dictionary of breeding values for the breeding program.
        gmod : dict
            Dictionary of genomic models for the breeding program.
        t_cur : int
            Current time in the breeding program.
        t_max : int
            Deadline time for the breeding program.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            A tuple of length 5: ``(genome, geno, pheno, bval, gmod)``

            Where:

            - ``genome`` is a ``dict`` of genomes for the breeding program.
            - ``geno`` is a ``dict`` of genotypes for the breeding program.
            - ``pheno`` is a ``dict`` of phenotypes for the breeding program.
            - ``bval`` is a ``dict`` of breeding values for the breeding program.
            - ``gmod`` is a ``dict`` of genomic models for the breeding program.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_SurvivorSelectionOperator(v: object, vname: str) -> None:
    """
    Check if object is of type SurvivorSelectionOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, SurvivorSelectionOperator):
        raise TypeError("'%s' must be a SurvivorSelectionOperator." % vname)
