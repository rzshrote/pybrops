"""
Module defining interfaces and associated error checking routines for mating
mating protocols.
"""

from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from typing import Optional
from typing import Union
import numpy

from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class MatingProtocol(
        metaclass = ABCMeta,
    ):
    """
    Abstract class for mating protocols.

    The purpose of this abstract class is to provide functionality for:
        1) Mating simulation and progeny generation from genotype matrices.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################
    @property
    @abstractmethod
    def nparent(self) -> Integral:
        """Number of parents the mating protocol requires."""
        raise NotImplementedError("property is abstract")
    @nparent.setter
    @abstractmethod
    def nparent(self, value: Integral) -> None:
        """Set number of parents the mating protocol requires."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################
    @abstractmethod
    def mate(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            xconfig: numpy.ndarray, 
            nmating: Union[Integral,numpy.ndarray], 
            nprogeny: Union[Integral,numpy.ndarray], 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> PhasedGenotypeMatrix:
        """
        Mate individuals according to a mating scheme.

        Parameters
        ----------
        pgmat : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of parental candidates.
        xconfig : numpy.ndarray
            Array of shape ``(ncross,nparent)`` containing indices specifying a cross
            configuration. Each index corresponds to an individual in ``pgmat``.
        nmating : Integral, numpy.ndarray
            Number of matings of the cross configuration per cross pattern.
            Relevant in situations with heterozygous parents.
        nprogeny : Integral, numpy.ndarray
            Number of progeny to generate per mating.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If ``dict``, write to dict (may overwrite previously defined fields).
            If ``None``, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : PhasedGenotypeMatrix
            A PhasedGenotypeMatrix of progeny.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_MatingProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type MatingProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, MatingProtocol):
        raise TypeError("'%s' must be a MatingProtocol." % vname)
