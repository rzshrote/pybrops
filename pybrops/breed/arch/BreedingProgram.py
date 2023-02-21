"""
Module defining interfaces and associated error checking routines for
germplasm bank representation.
"""

from typing import Any
from pybrops.breed.arch.BreedingNode import BreedingNode

class BreedingProgram(BreedingNode):
    """
    Abstract class defining a breeding program.

    The purpose of this abstract class is to define functionality for:
        1) Container storage for the breeding program.
        2) Contain breeding operators used in the breeding program.
        3) Initialization routines for the breeding program.
        4) Advancement of simulations for the breeding program.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class BreedingProgram.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingProgram, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############ Starting condition containers #############
    @property
    def start_genome(self) -> Any:
        """Starting genomes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_genome.setter
    def start_genome(self, value: Any) -> None:
        """Set starting genomes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")
    @start_genome.deleter
    def start_genome(self) -> None:
        """Delete starting genomes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def start_geno(self) -> Any:
        """Starting genotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_geno.setter
    def start_geno(self, value: Any) -> None:
        """Set starting genotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")
    @start_geno.deleter
    def start_geno(self) -> None:
        """Delete starting genotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def start_pheno(self) -> Any:
        """Starting phenotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_pheno.setter
    def start_pheno(self, value: Any) -> None:
        """Set starting phenotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")
    @start_pheno.deleter
    def start_pheno(self) -> None:
        """Delete starting phenotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def start_bval(self) -> Any:
        """Starting breeding values for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_bval.setter
    def start_bval(self, value: Any) -> None:
        """Set starting breeding values for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")
    @start_bval.deleter
    def start_bval(self) -> None:
        """Delete starting breeding values for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def start_gmod(self) -> Any:
        """Starting genomic models for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_gmod.setter
    def start_gmod(self, value: Any) -> None:
        """Set starting genomic models for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")
    @start_gmod.deleter
    def start_gmod(self) -> None:
        """Delete starting genomic models for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    ######### Breeding program operator properties #########
    @property
    def initop(self) -> Any:
        """Initialization operator."""
        raise NotImplementedError("property is abstract")
    @initop.setter
    def initop(self, value: Any) -> None:
        """Set the initialization operator"""
        raise NotImplementedError("property is abstract")
    @initop.deleter
    def initop(self) -> None:
        """Delete the initialization operator"""
        raise NotImplementedError("property is abstract")

    @property
    def pselop(self) -> Any:
        """Parent selection operator."""
        raise NotImplementedError("property is abstract")
    @pselop.setter
    def pselop(self, value: Any) -> None:
        """Set the parent selection operator"""
        raise NotImplementedError("property is abstract")
    @pselop.deleter
    def pselop(self) -> None:
        """Delete the parent selection operator"""
        raise NotImplementedError("property is abstract")

    @property
    def mateop(self) -> Any:
        """Mating operator."""
        raise NotImplementedError("property is abstract")
    @mateop.setter
    def mateop(self, value: Any) -> None:
        """Set the mating operator"""
        raise NotImplementedError("property is abstract")
    @mateop.deleter
    def mateop(self) -> None:
        """Delete the mating operator"""
        raise NotImplementedError("property is abstract")

    @property
    def evalop(self) -> Any:
        """Evaluation operator."""
        raise NotImplementedError("property is abstract")
    @evalop.setter
    def evalop(self, value: Any) -> None:
        """Set the evaluation operator"""
        raise NotImplementedError("property is abstract")
    @evalop.deleter
    def evalop(self) -> None:
        """Delete the evaluation operator"""
        raise NotImplementedError("property is abstract")

    @property
    def sselop(self) -> Any:
        """Survivor selection operator."""
        raise NotImplementedError("property is abstract")
    @sselop.setter
    def sselop(self, value: Any) -> None:
        """Set the survivor selection operator"""
        raise NotImplementedError("property is abstract")
    @sselop.deleter
    def sselop(self) -> None:
        """Delete the survivor selection operator"""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############# Initialize breeding program ##############
    def initialize(self, **kwargs: dict):
        """
        Initialize the breeding program with genotypes, phenotypes, and genomic
        models.
        """
        raise NotImplementedError("method is abstract")

    def is_initialized(self, **kwargs: dict):
        """
        Return whether or not the BreedingProgram has been initialized with a
        starting set of conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : boolean
            True if the BreedingProgram has been initialized.
            False if the BreedingProgram has not been initialized.
        """
        raise NotImplementedError("method is abstract")

    ############# Population evolution methods #############
    def reset(self, **kwargs: dict):
        """
        Reset the evolution of the breeding program back to starting conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def advance(self, ngen, lbook, **kwargs: dict):
        """
        Advance the breeding program by a specified number of generations.

        Parameters
        ----------
        ngen : int
            Number of generations to advance the BreedingProgram.
        lbook : Logbook
            Logbook into which to write statistics.
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def evolve(self, nrep, ngen, lbook, **kwargs: dict):
        """
        Evolve the breeding program for a set number of replications and
        generations. The BreedingProgram is restarted using the starting geno,
        bval, gmod containers.

        Parameters
        ----------
        nrep : int
            Number of evolution replicates.
        ngen : int
            Number of generations to evolve the population for each replicate.
            Note that this does not modify 't_max'.
        lbook : Logbook
            Logbook into which to write statistics.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingProgram(v: Any) -> bool:
    """
    Determine whether an object is a BreedingProgram.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingProgram object instance.
    """
    return isinstance(v, BreedingProgram)

def check_is_BreedingProgram(v: Any, varname: str) -> None:
    """
    Check if object is of type BreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingProgram):
        raise TypeError("'%s' must be a BreedingProgram." % varname)
