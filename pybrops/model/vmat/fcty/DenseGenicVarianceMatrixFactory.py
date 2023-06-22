from typing import Any
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.DenseGenicVarianceMatrix import DenseGenicVarianceMatrix
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix
from pybrops.model.vmat.fcty.GenicVarianceMatrixFactory import GenicVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DenseGenicVarianceMatrixFactory(GenicVarianceMatrixFactory):
    """
    docstring for DenseGenicVarianceMatrixFactory.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, **kwargs: dict):
        """
        Constructor for DenseGenicVarianceMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseGenicVarianceMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmod(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            nprogeny: int, 
            **kwargs: dict
        ) -> DenseGenicVarianceMatrix:
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
        return DenseGenicVarianceMatrix.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            nprogeny = nprogeny, 
            **kwargs
        )



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
