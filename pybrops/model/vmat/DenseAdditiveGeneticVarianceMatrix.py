"""
Module implementing classes and associated error checking routines for matrices
storing dense additive genetic variance estimates.
"""

from typing import Any, Optional

import numpy
from pybrops.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
from pybrops.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix

# TODO: implement me
class DenseAdditiveGeneticVarianceMatrix(DenseGeneticVarianceMatrix,AdditiveGeneticVarianceMatrix):
    """
    A semi-concrete class for dense additive genetic variance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense additive genetic variance matrices.

    Methods responsible for estimating genetic variances from genomic models
    remain abstract and must be implemented by inheriting classes.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat: numpy.ndarray, taxa: Optional[numpy.ndarray] = None, taxa_grp: Optional[numpy.ndarray] = None, **kwargs: dict):
        """
        Constructor for the concrete class DenseAdditiveGeneticVarianceMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            Array used to construct the object.
        taxa : numpy.ndarray
            Taxa names.
        taxa_grp : numpy.ndarray
            Taxa groupings.
        kwargs : dict
            Additional keyword arguments.
        """
        super(DenseAdditiveGeneticVarianceMatrix, self).__init__(
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
def is_DenseAdditiveGeneticVarianceMatrix(v: Any) -> bool:
    """
    Determine whether an object is a ``DenseAdditiveGeneticVarianceMatrix``.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``DenseAdditiveGeneticVarianceMatrix`` object instance.
    """
    return isinstance(v, DenseAdditiveGeneticVarianceMatrix)

def check_is_DenseAdditiveGeneticVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``DenseAdditiveGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseAdditiveGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a DenseAdditiveGeneticVarianceMatrix".format(vname))
