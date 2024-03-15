"""
Module defining basal interfaces and associated error checking routines for
breeding nodes. Breeding nodes compose complex breeding programs in a graph-like
structure. They are points were germplasm and information are located.
"""


from abc import ABCMeta, abstractmethod


class BreedingNode(
        metaclass = ABCMeta,
    ):
    """
    Abstract class defining a breeding node. Breeding nodes compose complex
    breeding programs in a graph-like structure. They are points were germplasm
    and information are located.

    The purpose of this abstract class is to define functionality for:
        1) Container storage for germplasm and information.
        2) Time related information.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############ Program information containers ############
    @property
    @abstractmethod
    def genome(self) -> object:
        """Genomes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @genome.setter
    @abstractmethod
    def genome(self, value: object) -> None:
        """Set genomes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def geno(self) -> object:
        """Genotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @geno.setter
    @abstractmethod
    def geno(self, value: object) -> None:
        """Set genotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def pheno(self) -> object:
        """Phenotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @pheno.setter
    @abstractmethod
    def pheno(self, value: object) -> None:
        """Set phenotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def bval(self) -> object:
        """Breeding values for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @bval.setter
    @abstractmethod
    def bval(self, value: object) -> None:
        """Set breeding values for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def gmod(self) -> object:
        """Genomic models for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @gmod.setter
    @abstractmethod
    def gmod(self, value: object) -> None:
        """Set genomic models for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    ############# Generation number properties #############
    @property
    @abstractmethod
    def t_cur(self) -> int:
        """Current time of the BreedingNode."""
        raise NotImplementedError("property is abstract")
    @t_cur.setter
    @abstractmethod
    def t_cur(self, value: int) -> None:
        """Set the current time for the BreedingNode"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def t_max(self) -> int:
        """Maximum time of the BreedingNode."""
        raise NotImplementedError("property is abstract")
    @t_max.setter
    @abstractmethod
    def t_max(self, value: int) -> None:
        """Set the maximum time for the BreedingNode"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################



################################## Utilities ###################################
def check_is_BreedingNode(v: object, vname: str) -> None:
    """
    Check if object is of type BreedingNode. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingNode):
        raise TypeError("'%s' must be a BreedingNode." % vname)
