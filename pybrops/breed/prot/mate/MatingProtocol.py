"""
Module defining interfaces and associated error checking routines for mating
mating protocols.
"""

from numbers import Integral
from typing import Optional, Union
import numpy

from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class MatingProtocol:
    """
    Abstract class for mating protocols.

    The purpose of this abstract class is to provide functionality for:
        1) Mating simulation and progeny generation from genotype matrices.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for abstract class MatingProtocol.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments.
        """
        super(MatingProtocol, self).__init__()

    ############################ Object Properties #############################
    @property
    def nparent(self) -> Integral:
        """Number of parents the mating protocol requires."""
        raise NotImplementedError("property is abstract")
    @nparent.setter
    def nparent(self, value: Integral) -> None:
        """Set number of parents the mating protocol requires."""
        raise NotImplementedError("property is abstract")
    @nparent.deleter
    def nparent(self) -> None:
        """Delete number of parents the mating protocol requires."""
        raise NotImplementedError("property is abstract")
    

    ############################## Object Methods ##############################
    def mate(
            self, 
            pgmat: PhasedGenotypeMatrix, 
            xconfig: numpy.ndarray, 
            ncross: Union[int,numpy.ndarray], 
            nprogeny: Union[int,numpy.ndarray], 
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
            Array of indices specifying a cross configuration. Each index corresponds
            to an individual in ``pgmat``.
        ncross : numpy.ndarray
            Number of crosses to perform per cross pattern.
        nprogeny : numpy.ndarray
            Number of progeny to generate per cross.
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_MatingProtocol(v: object) -> bool:
    """
    Determine whether an object is a MatingProtocol.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a MatingProtocol object instance.
    """
    return isinstance(v, MatingProtocol)

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
