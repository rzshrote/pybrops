"""
Module defining interfaces and associated error checking routines for matrices
storing genic covariance estimates.
"""

__all__ = [
    "ProgenyGenicCovarianceMatrix",
    "check_is_ProgenyGenicCovarianceMatrix",
]

from abc import ABCMeta, abstractmethod
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.SquareTaxaSquareTraitMatrix import SquareTaxaSquareTraitMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class ProgenyGenicCovarianceMatrix(
        SquareTaxaSquareTraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for additive genetic covariance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of genic covariance from a genomic model.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ######## Expected parental genome contributions ########
    @property
    @abstractmethod
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring."""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ############################## Object Methods ##############################
    @classmethod
    @abstractmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int,
            **kwargs: dict
        ):
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : ProgenyGenicCovarianceMatrix
            A matrix of genic covariance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_ProgenyGenicCovarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``ProgenyGenicCovarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, ProgenyGenicCovarianceMatrix):
        raise TypeError("'{0}' must be a ProgenyGenicCovarianceMatrix".format(vname))
