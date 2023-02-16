"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense additive genic variance estimates.
"""

from typing import Any, Optional
import numpy
from pybrops.model.vmat.DenseGenicVarianceMatrix import DenseGenicVarianceMatrix
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix

# TODO: implement me
class DenseAdditiveGenicVarianceMatrix(DenseGenicVarianceMatrix,AdditiveGenicVarianceMatrix):
    """
    UNDER CONSTRUCTION!

    A semi-concrete class for dense additive genetic variance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense additive genetic variance matrices.

    Methods responsible for estimating genetic variances from genomic models
    remain abstract and must be implemented by inheriting classes.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ):
        super(DenseAdditiveGenicVarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    # from_gmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete

    # from_algmod
    # this method should remain abstract; it depends on the cross structure
    # maybe in the future if generic cross structures are implemented,
    # then this method can become concrete



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_DenseAdditiveGenicVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``DenseAdditiveGenicVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseAdditiveGenicVarianceMatrix):
        raise TypeError("'{0}' must be a DenseAdditiveGenicVarianceMatrix".format(vname))
