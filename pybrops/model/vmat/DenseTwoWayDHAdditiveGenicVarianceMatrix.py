"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense additive genic variance estimates calculated using two-way DH
formulae.
"""

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

    def __init__(self, arg):
        super(DenseTwoWayDHAdditiveGenicVarianceMatrix, self).__init__()
        self.arg = arg
