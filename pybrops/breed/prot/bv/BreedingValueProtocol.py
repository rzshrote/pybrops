"""
Module defining interfaces and associated error checking methods for breeding
value calculation protocols.
"""

from abc import ABCMeta
from abc import abstractmethod
from typing import Optional

import pandas

from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class BreedingValueProtocol(
        metaclass = ABCMeta,
    ):
    """
    Abstract class defining interfaces for breeding value calculation protocols.

    The purpose of this abstract class is to provide functionality for:
        1) Estimation of breeding values from phenotype values.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def estimate(
            self, 
            ptobj: pandas.DataFrame, 
            gtobj: GenotypeMatrix, 
            miscout: Optional[dict], 
            **kwargs: dict
        ) -> BreedingValueMatrix:
        """
        Estimate breeding values.

        Parameters
        ----------
        ptobj : pandas.DataFrame
            An object containing phenotype data. Must be a phenotype data frame.
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a genotype matrix.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            A matrix of breeding values.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_BreedingValueProtocol(v: object, vname: str) -> None:
    """
    Check if object is of type BreedingValueProtocol. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingValueProtocol):
        raise TypeError("'%s' must be a BreedingValueProtocol." % vname)
