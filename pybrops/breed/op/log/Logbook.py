"""
Module defining interfaces and associated error checking routines for
breeding program logbook operators.
"""

__all__ = [
    "Logbook",
    "check_is_Logbook"
]

from typing import Any

class Logbook:
    """
    Abstract class defining interfaces for logging statistics about an entire
    breeding program.

    The purpose of this abstract class is to provide functionality for:
        1) Logging data on an entire breeding program.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class Logbook.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(Logbook, self).__init__()

    ############################ Object Properties #############################
    @property
    def data(self) -> Any:
        """Logbook data."""
        raise NotImplementedError("property is abstract")
    @data.setter
    def data(self, value: Any) -> None:
        """Set logbook data."""
        raise NotImplementedError("property is abstract")

    @property
    def rep(self) -> int:
        """Replicate number."""
        raise NotImplementedError("property is abstract")
    @rep.setter
    def rep(self, value: int) -> None:
        """Set replicate number."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################
    def log_initialize(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs: dict):
        """
        Record information directly after 'InitializationOperator.initialize'
        is called.

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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_pselect(self, mcfg, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs: dict):
        """
        Record information directly after 'ParentSelectionOperator.pselect'
        is called.

        Parameters
        ----------
        mcfg : dict
            Dictionary of mating configurations for the breeding program.
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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_mate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs: dict):
        """
        Record information directly after 'MatingOperator.mate' is called.

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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_evaluate(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs: dict):
        """
        Record information directly after 'EvaluationOperator.evaluate' is
        called.

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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def log_sselect(self, genome, geno, pheno, bval, gmod, t_cur, t_max, **kwargs: dict):
        """
        Record information directly after 'SurvivorSelectionOperator.sselect'
        is called.

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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def reset(self):
        """
        Reset Logbook internals.
        """
        raise NotImplementedError("method is abstract")

    def write(self, fname):
        """
        Write Logbook to file

        Parameters
        ----------
        fname : str
            File name to which to write file.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_Logbook(v: object, vname: str) -> None:
    """
    Check if object is of type Logbook. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, Logbook):
        raise TypeError("'%s' must be a Logbook." % vname)
