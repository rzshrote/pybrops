from abc import ABCMeta, abstractmethod
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class AdditiveGenicVarianceMatrixFactory(GenicVarianceMatrixFactory,metaclass=ABCMeta):
    """
    Abstract factory class for producing AdditiveGenicVarianceMatrix objects.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def from_algmod(
            self, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            mem: int,
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
        out : GeneticVarianceMatrix
            A matrix of additive genic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_AdditiveGenicVarianceMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type ``AdditiveGenicVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveGenicVarianceMatrixFactory):
        raise TypeError("'{0}' must be a AdditiveGenicVarianceMatrixFactory".format(vname))
