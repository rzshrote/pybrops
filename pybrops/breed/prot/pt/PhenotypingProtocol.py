"""
Module defining interfaces and associated protocols for phenotyping protocols.
"""

import numbers
from typing import Any, Union

import numpy

from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame


class PhenotypingProtocol:
    """
    Abstract class defining interfaces for phenotyping protocols.

    The purpose of this abstract class is to provide functionality for:
        1) Genomic model metadata.
        2) Phenotype simulation.
        3) Manipulation and setting of environmental variance metadata.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class PhenotypingProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(PhenotypingProtocol, self).__init__()

    ############################ Object Properties #############################

    ############### Genomic Model Properties ###############
    @property
    def gpmod(self) -> GenomicModel:
        """Genomic prediction model."""
        raise NotImplementedError("property is abstract")
    @gpmod.setter
    def gpmod(self, value: GenomicModel) -> None:
        """Set genomic prediction model"""
        raise NotImplementedError("property is abstract")
    @gpmod.deleter
    def gpmod(self) -> None:
        """Delete genomic prediction model"""
        raise NotImplementedError("property is abstract")

    ################ Stochastic Parameters #################
    @property
    def var_err(self) -> Any:
        """Error variance for each trait."""
        raise NotImplementedError("property is abstract")
    @var_err.setter
    def var_err(self, value: Any) -> None:
        """Set error variance"""
        raise NotImplementedError("property is abstract")
    @var_err.deleter
    def var_err(self) -> None:
        """Delete error variance"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################
    def phenotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: dict, 
            **kwargs: dict
        ) -> PhenotypeDataFrame:
        """
        Phenotype a set of genotypes using a genomic prediction model.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            Genomes of the individuals to phenotype.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhenotypeDataFrame
            A PhenotypeDataFrame containing phenotypes for individuals.
        """
        raise NotImplementedError("method is abstract")

    def set_h2(
            self, 
            h2: Union[numbers.Number,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")

    def set_H2(
            self, 
            H2: Union[numbers.Number,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : float, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypingProtocol(v: object) -> bool:
    """
    Determine whether an object is a PhenotypingProtocol.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a PhenotypingProtocol object instance.
    """
    return isinstance(v, PhenotypingProtocol)

def check_is_PhenotypingProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type PhenotypingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, PhenotypingProtocol):
        raise TypeError("'%s' must be a PhenotypingProtocol." % vname)
