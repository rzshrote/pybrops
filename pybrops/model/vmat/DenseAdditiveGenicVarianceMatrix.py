from pybrops.model.vmat.DenseGenicVarianceMatrix import DenseGenicVarianceMatrix
from pybrops.model.vmat.AdditiveGenicVarianceMatrix import AdditiveGenicVarianceMatrix

# TODO: implement me
class DenseAdditiveGenicVarianceMatrix(DenseGenicVarianceMatrix,AdditiveGenicVarianceMatrix):
    """docstring for DenseAdditiveGenicVarianceMatrix."""

    def __init__(self, arg):
        super(DenseAdditiveGenicVarianceMatrix, self).__init__()
        self.arg = arg
