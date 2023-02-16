from typing import Any
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix
from pybrops.model.vmat.GeneticVarianceMatrixFactory import GeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix


class AdditiveGeneticVarianceMatrixFactory(GeneticVarianceMatrixFactory):
    """
    Abstract factory class for producing AdditiveGeneticVarianceMatrix objects.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for AdditiveGeneticVarianceMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(AdditiveGeneticVarianceMatrixFactory, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def from_algmod(
            self, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
            nprogeny: int, 
            nself: int, 
            gmapfn: GeneticMapFunction, 
            mem: int,
            **kwargs: dict
        ) -> AdditiveGeneticVarianceMatrix:
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
        out : GeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_AdditiveGeneticVarianceMatrixFactory(v: Any, vname: str) -> None:
    """
    Check if object is of type ``AdditiveGeneticVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveGeneticVarianceMatrixFactory):
        raise TypeError("'{0}' must be a AdditiveGeneticVarianceMatrixFactory".format(vname))
