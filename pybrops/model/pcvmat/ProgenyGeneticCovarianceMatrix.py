"""
Module defining interfaces and error checking routines for matrices storing progeny genetic covariance-covariance estimates.
"""

from abc import ABCMeta, abstractmethod
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.SquareTaxaSquareTraitMatrix import SquareTaxaSquareTraitMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class ProgenyGeneticCovarianceMatrix(
        SquareTaxaSquareTraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for matrices storing progeny covariances for specific crosses.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of genetic covariances from a genomic model.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ######## Expected parental genome contributions ########
    @property
    @abstractmethod
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring from each parent."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> 'ProgenyGeneticCovarianceMatrix':
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
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
            Genetic map function with which to calculate recombination probabilities.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : ProgenyGeneticCovarianceMatrix
            A matrix of genetic covariance estimations.
        """
        raise NotImplementedError("method is abstract")

    ############################## Static Methods ##############################



################################## Utilities ###################################
def check_is_ProgenyGeneticCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``ProgenyGeneticCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ProgenyGeneticCovarianceMatrix):
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,ProgenyGeneticCovarianceMatrix.__name__,type(v).__name__))
