"""
Module defining interfaces and associated error checking routines for
germplasm bank representation.
"""

__all__ = [
    "BreedingProgram",
    "check_is_BreedingProgram"
]

from abc import ABCMeta, abstractmethod
from pybrops.breed.arch.BreedingNode import BreedingNode
from pybrops.breed.op.log.Logbook import Logbook

class BreedingProgram(BreedingNode,metaclass=ABCMeta):
    """
    Abstract class defining a breeding program.

    The purpose of this abstract class is to define functionality for:
        1) Container storage for the breeding program.
        2) Contain breeding operators used in the breeding program.
        3) Initialization routines for the breeding program.
        4) Advancement of simulations for the breeding program.
    """

    ########################## Special Object Methods ##########################
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

    ############################ Object Properties #############################

    ############ Starting condition containers #############
    @property
    @abstractmethod
    def start_genome(self) -> object:
        """Starting genomes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_genome.setter
    @abstractmethod
    def start_genome(self, value: object) -> None:
        """Set starting genomes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def start_geno(self) -> object:
        """Starting genotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_geno.setter
    @abstractmethod
    def start_geno(self, value: object) -> None:
        """Set starting genotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def start_pheno(self) -> object:
        """Starting phenotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_pheno.setter
    @abstractmethod
    def start_pheno(self, value: object) -> None:
        """Set starting phenotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def start_bval(self) -> object:
        """Starting breeding values for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_bval.setter
    @abstractmethod
    def start_bval(self, value: object) -> None:
        """Set starting breeding values for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def start_gmod(self) -> object:
        """Starting genomic models for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @start_gmod.setter
    @abstractmethod
    def start_gmod(self, value: object) -> None:
        """Set starting genomic models for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    ######### Breeding program operator properties #########
    @property
    @abstractmethod
    def initop(self) -> object:
        """Initialization operator."""
        raise NotImplementedError("property is abstract")
    @initop.setter
    @abstractmethod
    def initop(self, value: object) -> None:
        """Set the initialization operator"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def pselop(self) -> object:
        """Parent selection operator."""
        raise NotImplementedError("property is abstract")
    @pselop.setter
    @abstractmethod
    def pselop(self, value: object) -> None:
        """Set the parent selection operator"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def mateop(self) -> object:
        """Mating operator."""
        raise NotImplementedError("property is abstract")
    @mateop.setter
    @abstractmethod
    def mateop(self, value: object) -> None:
        """Set the mating operator"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def evalop(self) -> object:
        """Evaluation operator."""
        raise NotImplementedError("property is abstract")
    @evalop.setter
    @abstractmethod
    def evalop(self, value: object) -> None:
        """Set the evaluation operator"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def sselop(self) -> object:
        """Survivor selection operator."""
        raise NotImplementedError("property is abstract")
    @sselop.setter
    @abstractmethod
    def sselop(self, value: object) -> None:
        """Set the survivor selection operator"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ############# Initialize breeding program ##############
    @abstractmethod
    def initialize(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Initialize the breeding program with genotypes, phenotypes, and genomic
        models.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def is_initialized(
            self, 
            **kwargs: dict
        ) -> bool:
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
    @abstractmethod
    def reset(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Reset the evolution of the breeding program back to starting conditions.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def advance(
            self, 
            ngen: int, 
            lbook: Logbook, 
            **kwargs: dict
        ) -> None:
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

    @abstractmethod
    def evolve(
            self, 
            nrep: int, 
            ngen: int, 
            lbook: Logbook, 
            **kwargs: dict
        ) -> None:
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



################################## Utilities ###################################
def check_is_BreedingProgram(v: object, vname: str) -> None:
    """
    Check if object is of type BreedingProgram. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingProgram):
        raise TypeError("'%s' must be a BreedingProgram." % vname)
