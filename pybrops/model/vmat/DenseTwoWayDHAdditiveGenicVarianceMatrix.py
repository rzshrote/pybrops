"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense additive genic variance estimates calculated using two-way DH
formulae.
"""

from typing import Optional
import numpy
from pybrops.model.vmat.DenseAdditiveGenicVarianceMatrix import DenseAdditiveGenicVarianceMatrix

# TODO: implement me
class DenseTwoWayDHAdditiveGenicVarianceMatrix(DenseAdditiveGenicVarianceMatrix):
    """
    UNDER CONSTRUCTION!
    A concrete class for dense additive genic variance matrices calculated
    for two-way DH progenies.

    The purpose of this concrete class is to implement functionality for:
        1) Genic variance estimation for two-way DH progenies.
    """

    def __init__(
            self,
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ):
        super(DenseTwoWayDHAdditiveGenicVarianceMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            **kwargs
        )
