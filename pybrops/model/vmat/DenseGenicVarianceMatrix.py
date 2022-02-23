"""
UNDER CONSTRUCTION!
Module implementing classes and associated error checking routines for matrices
storing dense genic variance estimates.
"""

from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix

# TODO: implement me
class DenseGenicVarianceMatrix(DenseSquareTaxaMatrix,GenicVarianceMatrix):
    """
    UNDER CONSTRUCTION!
    A semi-concrete class for dense genic variance matrices.

    The purpose of this semi-concrete class is to provide functionality for:
        1) Object construction of dense genic variance matrices.

    Methods responsible for estimating genic variances from genomic models
    remain abstract and must be implemented by inheriting classes.
    """

    def __init__(self, arg):
        super(DenseGenicVarianceMatrix, self).__init__()
        self.arg = arg
