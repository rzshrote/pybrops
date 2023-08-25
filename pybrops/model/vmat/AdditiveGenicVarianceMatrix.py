"""
Module defining interfaces and associated error checking routines for matrices
storing additive genic variance estimates.
"""

from abc import ABCMeta, abstractmethod
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class AdditiveGenicVarianceMatrix(GenicVarianceMatrix,metaclass=ABCMeta):
    """
    An abstract class for additive genetic variance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of additive genic variance from an additive linear
           genomic model.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############################## Object Methods ##############################
    @classmethod
    @abstractmethod
    def from_algmod(
            cls, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            mem: int
        ) -> 'AdditiveGenicVarianceMatrix':
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : GenicVarianceMatrix
            A matrix of additive genic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_AdditiveGenicVarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``AdditiveGenicVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveGenicVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGenicVarianceMatrix".format(vname))
