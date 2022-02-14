from pybrops.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybrops.model.vmat.GenicVarianceMatrix import GenicVarianceMatrix

# TODO: implement me
class DenseGenicVarianceMatrix(DenseSquareTaxaMatrix,GenicVarianceMatrix):
    """docstring for DenseGenicVarianceMatrix."""

    def __init__(self, arg):
        super(DenseGenicVarianceMatrix, self).__init__()
        self.arg = arg
