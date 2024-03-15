"""
Module defining interfaces and associated protocols for phenotyping protocols.
"""

__all__ = [
    "PhenotypingProtocol",
    "check_is_PhenotypingProtocol",
]

from abc import ABCMeta, abstractmethod
from numbers import Real
from typing import Union
import numpy
import pandas
from pybrops.core.io.Copyable import Copyable
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class PhenotypingProtocol(
        Copyable,
        HDF5InputOutput,
        metaclass=ABCMeta
    ):
    """
    Abstract class defining interfaces for phenotyping protocols.

    The purpose of this abstract class is to provide functionality for:
        1) Genomic model metadata.
        2) Phenotype simulation.
        3) Manipulation and setting of environmental variance metadata.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############### Genomic Model Properties ###############
    @property
    @abstractmethod
    def gpmod(self) -> GenomicModel:
        """Genomic prediction model."""
        raise NotImplementedError("property is abstract")
    @gpmod.setter
    @abstractmethod
    def gpmod(self, value: GenomicModel) -> None:
        """Set genomic prediction model"""
        raise NotImplementedError("property is abstract")

    ################ Stochastic Parameters #################
    @property
    @abstractmethod
    def var_err(self) -> object:
        """Pure error variance for each trait"""
        raise NotImplementedError("property is abstract")
    @var_err.setter
    @abstractmethod
    def var_err(self, value: object) -> None:
        """Set the pure error variance for each trait"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################
    @abstractmethod
    def phenotype(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            miscout: dict, 
            **kwargs: dict
        ) -> pandas.DataFrame:
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
        out : pandas.DataFrame
            A pandas.DataFrame containing phenotypes for individuals.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def set_h2(
            self, 
            h2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        h2 : Real, numpy.ndarray
            Narrow sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def set_H2(
            self, 
            H2: Union[Real,numpy.ndarray], 
            pgmat: PhasedGenotypeMatrix, 
            **kwargs: dict
        ) -> None:
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        H2 : Real, numpy.ndarray
            Broad sense heritability.
        pgmat : PhasedGenotypeMatrix
            Founder genotypes.
        kwargs : dict
            Additional keyword arguments
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
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
