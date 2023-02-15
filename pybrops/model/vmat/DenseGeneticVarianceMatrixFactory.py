from typing import Any
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class DenseGeneticVarianceMatrixFactory(GeneticVarianceMatrixFactory):
    """
    docstring for DenseGeneticVarianceMatrixFactory.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        """
        Constructor for DenseGeneticVarianceMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseGeneticVarianceMatrixFactory, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def from_gmod(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> GeneticVarianceMatrix:
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
        out : GeneticVarianceMatrix
            A matrix of genetic variance estimations.
        """
        return DenseGeneticVarianceMatrix.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            ncross = ncross, 
            nprogeny = nprogeny, 
            nself = nself, 
            gmapfn = gmapfn, 
            **kwargs
        )



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_GeneticVarianceMatrixFactory(v: Any, vname: str) -> None:
    """
    Check if object is of type ``GeneticVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, GeneticVarianceMatrixFactory):
        raise TypeError("'{0}' must be a GeneticVarianceMatrixFactory".format(vname))
