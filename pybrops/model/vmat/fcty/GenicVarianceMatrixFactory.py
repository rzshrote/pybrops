from abc import ABCMeta
from abc import abstractmethod
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class GenicVarianceMatrixFactory(
        metaclass = ABCMeta,
    ):
    """
    Abstract class for GenicVarianceMatrix factory classes.

    The purpose of this abstract interface is to provide functionality for:
        1) Construction of genetic variance matrices.
    """

    ########################## Special Object Methods ##########################

    ############################## Object Methods ##############################
    @abstractmethod
    def from_gmod(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            **kwargs: dict
        ) -> GenicVarianceMatrix:
        """
        Estimate genetic variances from a GenomicModel and PhasedGenotypeMatrix.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : int
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
        out : GenicVarianceMatrix
            A matrix of genetic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GenicVarianceMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type ``GenicVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, GenicVarianceMatrixFactory):
        raise TypeError("'{0}' must be a GenicVarianceMatrixFactory".format(vname))
