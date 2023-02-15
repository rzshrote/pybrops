"""
Module defining interfaces and associated error checking routines for matrices
storing genetic variance estimates.
"""

from typing import Any
from pybrops.core.mat.SquareTaxaMatrix import SquareTaxaMatrix
from pybrops.model.gmod.GenomicModel import GenomicModel
from pybrops.popgen.gmap.GeneticMapFunction import GeneticMapFunction
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class GeneticVarianceMatrix(SquareTaxaMatrix):
    """
    An abstract class for additive genetic variance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of genetic variance from a genomic model.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class GeneticVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(GeneticVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ######## Expected parental genome contributions ########
    @property
    def epgc(self) -> tuple:
        """Expected parental genome contribution to the offspring."""
        raise NotImplementedError("property is abstract")
    @epgc.getter
    def epgc(self) -> tuple:
        """Get a tuple of the expected parental genome contributions."""
        raise NotImplementedError("property is abstract")
    @epgc.setter
    def epgc(self, value: tuple) -> None:
        """Set a tuple of the expected parental genome contributions."""
        raise NotImplementedError("property is abstract")    
    @epgc.deleter
    def epgc(self) -> None:
        """Delete the expected parental genome contributions tuple."""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def from_gmod(
            cls, 
            gmod: GenomicModel, 
            pgmat: PhasedGenotypeMatrix, 
            ncross: int, 
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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GeneticVarianceMatrix(v: Any) -> bool:
    """
    Determine whether an object is a ``GeneticVarianceMatrix``.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``GeneticVarianceMatrix`` object instance.
    """
    return isinstance(v, GeneticVarianceMatrix)

def check_is_GeneticVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``GeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, GeneticVarianceMatrix):
        raise TypeError("'{0}' must be a GeneticVarianceMatrix".format(vname))
