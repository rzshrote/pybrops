"""
Module defining interfaces and associated error checking routines for matrices
storing genic variance estimates.
"""

__all__ = [
    "GenicVarianceMatrix",
    "check_is_GenicVarianceMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class GenicVarianceMatrix(
        SquareTaxaMatrix,
        TraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for additive genetic variance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of genic variance from a genomic model.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ######## Expected parental genome contributions ########
    @property
    @abstractmethod
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring."""
        raise NotImplementedError("property is abstract")
    @epgc.setter
    @abstractmethod
    def epgc(self, value: tuple) -> None:
        """Set a tuple of the expected parental genome contributions."""
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
        ) -> 'GenicVarianceMatrix':
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GenicVarianceMatrix
            A matrix of genic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GenicVarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``GenicVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, GenicVarianceMatrix):
        raise TypeError("'{0}' must be a GenicVarianceMatrix".format(vname))
