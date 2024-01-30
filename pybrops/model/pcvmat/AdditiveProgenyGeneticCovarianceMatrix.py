"""
Module defining interfaces and associated error checking routines for matrices
storing additive genetic covariance estimates.
"""

from abc import ABCMeta, abstractmethod
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.pcvmat.ProgenyGeneticCovarianceMatrix import ProgenyGeneticCovarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class AdditiveProgenyGeneticCovarianceMatrix(
        ProgenyGeneticCovarianceMatrix,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for additive genetic covariance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of additive genetic covariance from an additive linear
           genomic model.
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
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            mem: int,
            **kwargs: dict
        ) -> 'AdditiveProgenyGeneticCovarianceMatrix':
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : int
            Number of cross patterns to simulate for genetic covariance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            covariance.
        nself : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : int
            Memory chunk size to use during matrix operations.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : ProgenyGeneticCovarianceMatrix
            A matrix of additive genetic covariance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_AdditiveProgenyGeneticCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``AdditiveProgenyGeneticCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveProgenyGeneticCovarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveProgenyGeneticCovarianceMatrix".format(vname))
