"""
Module defining interfaces and associated error checking routines for matrices
storing additive genetic variance estimates.
"""

from typing import Any
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix

class AdditiveGeneticVarianceMatrix(GeneticVarianceMatrix):
    """
    An abstract class for additive genetic variance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of additive genetic variance from an additive linear
           genomic model.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class AdditiveGeneticVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(AdditiveGeneticVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    @classmethod
    def from_algmod(cls, algmod, pgmat, ncross, nprogeny, s, gmapfn, mem):
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
        s : int
            Number of selfing generations post-cross pattern before 'nprogeny'
            individuals are simulated.
        gmapfn : GeneticMapFunction
            GeneticMapFunction to use to estimate covariance induced by
            recombination.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : GeneticVarianceMatrix
            A matrix of additive genetic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_AdditiveGeneticVarianceMatrix(v: Any) -> bool:
    """
    Determine whether an object is a ``AdditiveGeneticVarianceMatrix``.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``v`` is a ``AdditiveGeneticVarianceMatrix`` object instance.
    """
    return isinstance(v, AdditiveGeneticVarianceMatrix)

def check_is_AdditiveGeneticVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``AdditiveGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGeneticVarianceMatrix".format(vname))
