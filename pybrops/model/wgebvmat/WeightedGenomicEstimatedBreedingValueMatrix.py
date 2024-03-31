"""
Module defining interfaces and error checking routines for weighted genomic 
estimated breeding value matrices.
"""

__all__ = [
    "WeightedGenomicEstimatedBreedingValueMatrix",
    "check_is_WeightedGenomicEstimatedBreedingValueMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class WeightedGenomicEstimatedBreedingValueMatrix(
        BreedingValueMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for storing weighted genomic estimated breeding value matrices.

    The purpose of this abstract interface is to provide functionality for:

        1) Estimation of wGEBVs from a genomic model.
    """
    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_algmod(
            cls,
            algmod: AdditiveLinearGenomicModel,
            gmat: GenotypeMatrix,
            **kwargs: dict
        ) -> 'WeightedGenomicEstimatedBreedingValueMatrix':
        """
        Estimate Weighted Genomic Estimated Breeding Values (wGEBVs) from an 
        ``AdditiveLinearGenomicModel``.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            An ``AdditiveLinearGenomicModel`` with which to estimate wGEBVs.
        
        gmat : GenotypeMatrix
            Genotypes for which to estimate wGEBVs.
        
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : WeightedGenomicEstimatedBreedingValueMatrix
            A ``WeightedGenomicEstimatedBreedingValueMatrix`` object containing 
            wGEBVs.
        """
        raise NotImplementedError("classmethod is abstract")
    
    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_WeightedGenomicEstimatedBreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if an object is of type ``WeightedGenomicEstimatedBreedingValueMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, WeightedGenomicEstimatedBreedingValueMatrix):
        raise TypeError("variable ``{0}`` must be of type ``{1}`` but received type ``{2}``".format(
                vname,
                WeightedGenomicEstimatedBreedingValueMatrix.__name__,
                type(v).__name__
            )
        )
