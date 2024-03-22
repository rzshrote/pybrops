"""
Module defining interfaces and associated error checking routines for matrices
storing additive genic covariance estimates.
"""

from abc import ABCMeta
from abc import abstractmethod
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.pcvmat.ProgenyGenicCovarianceMatrix import ProgenyGenicCovarianceMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class AdditiveProgenyGenicCovarianceMatrix(
        ProgenyGenicCovarianceMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for additive genetic covariance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of additive genic covariance from an additive linear
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
        ) -> 'AdditiveProgenyGenicCovarianceMatrix':
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
            covariance.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : ProgenyGenicCovarianceMatrix
            A matrix of additive genic covariance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_AdditiveProgenyGenicCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``AdditiveProgenyGenicCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveProgenyGenicCovarianceMatrix):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,AdditiveProgenyGenicCovarianceMatrix.__name__,type(v).__name__))
