"""
Module defining interfaces and error checking routines for expected maximum 
breeding value matrices.
"""

__all__ = [
    "ExpectedMaximumBreedingValueMatrix",
    "check_is_ExpectedMaximumBreedingValueMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from numbers import Integral
from typing import Union
import numpy
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class ExpectedMaximumBreedingValueMatrix(
        BreedingValueMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for storing expected maximum breeding value matrices.

    The purpose of this abstract interface is to provide functionality for:

        1) Estimation of EMBVs from a genomic model.
    """
    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_gmod(
            cls,
            gmod: GenomicModel,
            pgmat: PhasedGenotypeMatrix,
            nprogeny: Union[Integral,numpy.ndarray],
            nrep: Union[Integral,numpy.ndarray],
            **kwargs: dict
        ) -> 'ExpectedMaximumBreedingValueMatrix':
        """
        Estimate Expected Maximum Breeding Values (EMBVs) from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            A ``GenomicModel`` with which to estimate EMBVs.

        pgmat : PhasedGenotypeMatrix
            Genomes for which to estimate EMBVs.

        nprogeny : Integral, numpy.ndarray
            Number of DH progenies to create per individual.
            If ``Integral``, each taxon produces ``nprogeny`` DH progenies.
            If ``numpy.ndarray`` of shape ``(ntaxa,)``, each taxon produces the
            number of progenies corresponding to elements in the array.

        nrep : Integral, numpy.ndarray
            Number of sample replicates for which to sample ``nprogeny`` DH progenies.
            If ``Integral``, each taxon produces ``nrep`` sample replicates.
            If ``numpy.ndarray`` of shape ``(ntaxa,)``, each taxon produces the
            number of sample replicates corresponding to elements in the array.

        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : ExpectedMaximumBreedingValueMatrix
            An ``ExpectedMaximumBreedingValueMatrix`` object containing EMBVs.
        """
        raise NotImplementedError("classmethod is abstract")

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_ExpectedMaximumBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``ExpectedMaximumBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ExpectedMaximumBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                ExpectedMaximumBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
