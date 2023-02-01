"""
Module defining interfaces and associated error checking routines for
breeding program initialization operators.
"""

from typing import Any


class InitializationOperator:
    """
    Abstract class defining interfaces for the evaluation of an entire breeding
    program.

    The purpose of this abstract class is to provide functionality for:
        1) Initialization of an entire breeding program.
    """

    def __init__(self, **kwargs: dict) -> None:
        """
        Constructor for the abstract class InitializationOperator.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(InitializationOperator, self).__init__()

    def initialize(self, miscout, **kwargs: dict):
        """
        Initialize a breeding program.

        Parameters
        ----------
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_InitializationOperator(v: Any):
    """
    Determine whether an object is a InitializationOperator.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a InitializationOperator object instance.
    """
    return isinstance(v, InitializationOperator)

def check_is_InitializationOperator(v: Any, varname: str):
    """
    Check if object is of type InitializationOperator. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, InitializationOperator):
        raise TypeError("'%s' must be a InitializationOperator." % varname)
