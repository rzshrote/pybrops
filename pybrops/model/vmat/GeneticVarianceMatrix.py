"""
Module defining interfaces and associated error checking routines for matrices
storing genetic variance estimates.
"""

from abc import ABCMeta, abstractmethod
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.core.mat.TraitMatrix import TraitMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class GeneticVarianceMatrix(
        SquareTaxaMatrix,
        TraitMatrix,
        PandasInputOutput,
        CSVInputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for additive genetic variance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of genetic variance from a genomic model.
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

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nmating: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> 'GeneticVarianceMatrix':
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        nmating : int
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.
        gmapfn : GeneticMapFunction
            Genetic map function with which to calculate recombination probabilities.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of genetic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GeneticVarianceMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``GeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, GeneticVarianceMatrix):
        raise TypeError("'{0}' must be a GeneticVarianceMatrix".format(vname))
