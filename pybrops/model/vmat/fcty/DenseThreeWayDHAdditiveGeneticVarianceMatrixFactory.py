from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import AdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.DenseThreeWayDHAdditiveGeneticVarianceMatrix import DenseThreeWayDHAdditiveGeneticVarianceMatrix
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory(
        AdditiveGeneticVarianceMatrixFactory,
    ):
    """
    docstring for DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        pass

    ############################## Object Methods ##############################
    def from_gmod(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> DenseThreeWayDHAdditiveGeneticVarianceMatrix:
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
        out : DenseThreeWayDHAdditiveGeneticVarianceMatrix
            A matrix of genetic variance estimations.
        """
        return DenseThreeWayDHAdditiveGeneticVarianceMatrix.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            nmating = ncross, 
            nprogeny = nprogeny, 
            nself = nself, 
            gmapfn = gmapfn, 
            **kwargs
        )

    def from_algmod(
            self, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            mem: int = 1024,
            **kwargs: dict
        ) -> DenseThreeWayDHAdditiveGeneticVarianceMatrix:
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
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
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : int
            Memory chunk size to use during matrix operations.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : DenseThreeWayDHAdditiveGeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        return DenseThreeWayDHAdditiveGeneticVarianceMatrix.from_algmod(
            algmod = algmod, 
            pgmat = pgmat, 
            nmating = ncross, 
            nprogeny = nprogeny, 
            nself = nself, 
            gmapfn = gmapfn, 
            mem = mem, 
            **kwargs
        )



################################## Utilities ###################################
def check_is_DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory):
        raise TypeError("'{0}' must be a DenseThreeWayDHAdditiveGeneticVarianceMatrixFactory".format(vname))
