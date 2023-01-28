"""
Module implementing classes and associated error checking routines for matrices
storing dense genetic variance estimates.
"""

from typing import Any
from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix

# TODO: implement me
class DenseGeneticVarianceMatrix(DenseSquareTaxaMatrix,GeneticVarianceMatrix):
    """
    A semi-concrete class for dense genetic variance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense genetic variance matrices.

    Methods responsible for estimating genetic variances from genomic models
    remain abstract and must be implemented by inheriting classes.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, **kwargs):
        """
        Constructor for the concrete class DenseGeneticVarianceMatrix.

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
        super(DenseGeneticVarianceMatrix, self).__init__(
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGeneticVarianceMatrix(v: Any) -> bool:
    """
    Determine whether an object is a ``DenseGeneticVarianceMatrix``.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        ``True`` or ``False`` for whether ``obj`` is a ``DenseGeneticVarianceMatrix`` object instance.
    """
    return isinstance(v, DenseGeneticVarianceMatrix)

def check_is_DenseGeneticVarianceMatrix(v: Any, vname: str) -> None:
    """
    Check if object is of type ``DenseGeneticVarianceMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, DenseGeneticVarianceMatrix):
        raise TypeError("'{0}' must be a DenseGeneticVarianceMatrix".format(vname))
