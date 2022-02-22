"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense additive genic variance estimates.
"""

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

    def __init__(self, arg):
        super(DenseAdditiveGenicVarianceMatrix, self).__init__()
        self.arg = arg
