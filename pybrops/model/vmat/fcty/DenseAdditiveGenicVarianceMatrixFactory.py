from typing import Any
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix
from pybrops.model.vmat.DenseAdditiveGenicVarianceMatrix import DenseAdditiveGenicVarianceMatrix
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DenseAdditiveGenicVarianceMatrixFactory(GenicVarianceMatrixFactory):
    """
    Abstract factory class for producing AdditiveGenicVarianceMatrix objects.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseAdditiveGenicVarianceMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseAdditiveGenicVarianceMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmod(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            **kwargs: dict
        ) -> DenseAdditiveGenicVarianceMatrix:
        """
        Estimate genic variances from a GenomicModel and PhasedGenotypeMatrix.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genic variances.
        ncross : int
            Number of cross patterns to simulate for genic variance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genic
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
        out : GenicVarianceMatrix
            A matrix of genic variance estimations.
        """
        return DenseAdditiveGenicVarianceMatrix.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            nprogeny = nprogeny, 
            **kwargs
        )

    def from_algmod(
            self, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            mem: int = 1024,
            **kwargs: dict
        ) -> AdditiveGenicVarianceMatrix:
        """
        Estimate genic variances from a GenomicModel.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genic variances.
        ncross : int
            Number of cross patterns to simulate for genic variance
            estimation.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genic
            variance.
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
        out : GenicVarianceMatrix
            A matrix of additive genic variance estimations.
        """
        return DenseAdditiveGenicVarianceMatrix.from_algmod(
            algmod = algmod, 
            pgmat = pgmat, 
            nprogeny = nprogeny, 
            mem = mem
            **kwargs
        )



################################## Utilities ###################################
def check_is_DenseAdditiveGenicVarianceMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseAdditiveGenicVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseAdditiveGenicVarianceMatrixFactory):
        raise TypeError("'{0}' must be a DenseAdditiveGenicVarianceMatrixFactory".format(vname))
