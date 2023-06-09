"""
Module defining basal interfaces and associated error checking routines for
breeding nodes. Breeding nodes compose complex breeding programs in a graph-like
structure. They are points were germplasm and information are located.
"""

from typing import Any

class BreedingNode:
    """
    Abstract class defining a breeding node. Breeding nodes compose complex
    breeding programs in a graph-like structure. They are points were germplasm
    and information are located.

    The purpose of this abstract class is to define functionality for:
        1) Container storage for germplasm and information.
        2) Time related information.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class BreedingNode.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(BreedingNode, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############ Program information containers ############
    @property
    def genome(self) -> Any:
        """Genomes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @genome.setter
    def genome(self, value: Any) -> None:
        """Set genomes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def geno(self) -> Any:
        """Genotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @geno.setter
    def geno(self, value: Any) -> None:
        """Set genotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def pheno(self) -> Any:
        """Phenotypes for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @pheno.setter
    def pheno(self, value: Any) -> None:
        """Set phenotypes for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def bval(self) -> Any:
        """Breeding values for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @bval.setter
    def bval(self, value: Any) -> None:
        """Set breeding values for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    @property
    def gmod(self) -> Any:
        """Genomic models for individuals in the breeding program."""
        raise NotImplementedError("property is abstract")
    @gmod.setter
    def gmod(self, value: Any) -> None:
        """Set genomic models for individuals in the breeding program"""
        raise NotImplementedError("property is abstract")

    ############# Generation number properties #############
    @property
    def t_cur(self) -> int:
        """Current time of the BreedingNode."""
        raise NotImplementedError("property is abstract")
    @t_cur.setter
    def t_cur(self, value: int) -> None:
        """Set the current time for the BreedingNode"""
        raise NotImplementedError("property is abstract")

    @property
    def t_max(self) -> int:
        """Maximum time of the BreedingNode."""
        raise NotImplementedError("property is abstract")
    @t_max.setter
    def t_max(self, value: int) -> None:
        """Set the maximum time for the BreedingNode"""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
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
