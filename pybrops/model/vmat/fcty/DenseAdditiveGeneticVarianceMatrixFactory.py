"""
Module containing code for creating dense additive genetic 
variance matrices using a factory class.
"""

from numbers import Integral
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.model.vmat.fcty.AdditiveGeneticVarianceMatrixFactory import AdditiveGeneticVarianceMatrixFactory
from pybrops.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix
from pybrops.model.vmat.fcty.DenseGeneticVarianceMatrixFactory import DenseGeneticVarianceMatrixFactory
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class DenseAdditiveGeneticVarianceMatrixFactory(DenseGeneticVarianceMatrixFactory,AdditiveGeneticVarianceMatrixFactory):
    """
    docstring for DenseAdditiveGeneticVarianceMatrixFactory.
    """

    ########################## Special Object Methods ##########################
    def __init__(self, **kwargs: dict):
        """
        Constructor for DenseAdditiveGeneticVarianceMatrixFactory.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(DenseAdditiveGeneticVarianceMatrixFactory, self).__init__(**kwargs)

    ############################## Object Methods ##############################
    def from_gmod(
            self, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral, 
            gmapfn: GeneticMapFunction, 
            **kwargs: dict
        ) -> DenseAdditiveGeneticVarianceMatrix:
        """
        Estimate genetic variances from a GenomicModel and PhasedGenotypeMatrix.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : Integral
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
        return DenseAdditiveGeneticVarianceMatrix.from_gmod(
            gmod = gmod, 
            pgmat = pgmat, 
            ncross = ncross, 
            nprogeny = nprogeny, 
            nself = nself, 
            gmapfn = gmapfn, 
            **kwargs
        )

    def from_algmod(
            self, 
            algmod: AdditiveLinearGenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: Integral, 
            nprogeny: Integral, 
            nself: Integral, 
            gmapfn: GeneticMapFunction, 
            mem: Integral = 1024,
            **kwargs: dict
        ) -> DenseAdditiveGeneticVarianceMatrix:
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        algmod : AdditiveLinearGenomicModel
            AdditiveLinearGenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        ncross : Integral
            Number of cross patterns to simulate for genetic variance
            estimation.
        nprogeny : Integral
            Number of progeny to simulate per cross to estimate genetic
            variance.
        nself : Integral
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : Integral
            Memory chunk size to use during matrix operations.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        return DenseAdditiveGeneticVarianceMatrix.from_algmod(
            algmod = algmod, 
            pgmat = pgmat, 
            ncross = ncross, 
            nprogeny = nprogeny, 
            nself = nself, 
            gmapfn = gmapfn, 
            mem = mem, 
            **kwargs
        )



################################## Utilities ###################################
def check_is_DenseAdditiveGeneticVarianceMatrixFactory(v: object, vname: str) -> None:
    """
    Check if object is of type ``DenseAdditiveGeneticVarianceMatrixFactory``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseAdditiveGeneticVarianceMatrixFactory):
        raise TypeError("'{0}' must be a DenseAdditiveGeneticVarianceMatrixFactory".format(vname))
