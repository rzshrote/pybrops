"""
Module containing abstract class definitions for selection configurations
"""

from abc import ABCMeta, abstractmethod
from numbers import Integral
from typing import Union

import numpy

from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class SelectionConfiguration(metaclass=ABCMeta):
    """
    docstring for SelectionConfiguration.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################
    @property
    @abstractmethod
    def ncross(self) -> Integral:
        """Number of cross configurations to consider."""
        raise NotImplementedError("method is abstract")
    @ncross.setter
    @abstractmethod
    def ncross(self, value: Integral) -> None:
        """Set number of cross configurations to consider."""
        raise NotImplementedError("method is abstract")
    
    @property
    @abstractmethod
    def nparent(self) -> Integral:
        """Number of parents in a given cross configuration."""
        raise NotImplementedError("method is abstract")
    @nparent.setter
    @abstractmethod
    def nparent(self, value: Integral) -> None:
        """Set number of parents in a given cross configuration."""
        raise NotImplementedError("method is abstract")
    
    @property
    @abstractmethod
    def nmating(self) -> numpy.ndarray:
        """Number of times an individual cross configuration is executed."""
        raise NotImplementedError("method is abstract")
    @nmating.setter
    @abstractmethod
    def nmating(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of times an individual cross configuration is executed."""
        raise NotImplementedError("method is abstract")
    
    @property
    @abstractmethod
    def nprogeny(self) -> numpy.ndarray:
        """Number of progeny to derive from a mating event."""
        raise NotImplementedError("method is abstract")
    @nprogeny.setter
    @abstractmethod
    def nprogeny(self, value: Union[Integral,numpy.ndarray]) -> None:
        """Set number of progeny to derive from a mating event."""
        raise NotImplementedError("method is abstract")

    @property
    @abstractmethod
    def pgmat(self) -> PhasedGenotypeMatrix:
        """Genome matrix for the parental candidates."""
        raise NotImplementedError("method is abstract")
    @pgmat.setter
    @abstractmethod
    def pgmat(self, value: PhasedGenotypeMatrix) -> None:
        """Set genome matrix for the parental candidates."""
        raise NotImplementedError("method is abstract")
    
    @property
    @abstractmethod
    def xconfig(self) -> numpy.ndarray:
        """xconfig."""
        raise NotImplementedError("method is abstract")
    @xconfig.setter
    @abstractmethod
    def xconfig(self, value: numpy.ndarray) -> None:
        """Set xconfig."""
        raise NotImplementedError("method is abstract")
    
    ############################## Object Methods ##############################

    ############################## Class Methods ###############################

    ############################## Static Methods ##############################
