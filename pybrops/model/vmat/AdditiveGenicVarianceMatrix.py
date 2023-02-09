"""
Module defining interfaces and associated error checking routines for matrices
storing additive genic variance estimates.
"""

from typing import Any
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix

class AdditiveGenicVarianceMatrix(GenicVarianceMatrix):
    """
    An abstract class for additive genetic variance matrices.

    The purpose of this abstract interface is to provide functionality for:
        1) Estimation of additive genic variance from an additive linear
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
        Constructor for the abstract class AdditiveGenicVarianceMatrix.

        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments. Used for cooperative inheritance.
            Dictionary passing unused arguments to the parent class constructor.
        """
        super(AdditiveGenicVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    @classmethod
    def from_algmod(cls, algmod, pgmat, nprogeny, mem):
        """
        Estimate genetic variances from a GenomicModel.

        Parameters
        ----------
        gmod : GenomicModel
            GenomicModel with which to estimate genetic variances.
        pgmat : PhasedGenotypeMatrix
            Input genomes to use to estimate genetic variances.
        nprogeny : int
            Number of progeny to simulate per cross to estimate genetic
            variance.
        mem : int
            Memory chunk size to use during matrix operations.

        Returns
        -------
        out : GenicVarianceMatrix
            A matrix of additive genic variance estimations.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_AdditiveGenicVarianceMatrix(v: Any) -> bool:
    """
    Determine whether an object is a ``AdditiveGenicVarianceMatrix``.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``AdditiveGenicVarianceMatrix`` object instance.
    """
    return isinstance(v, AdditiveGenicVarianceMatrix)

def check_is_AdditiveGenicVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``AdditiveGenicVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, AdditiveGenicVarianceMatrix):
        raise TypeError("'{0}' must be a AdditiveGenicVarianceMatrix".format(vname))
