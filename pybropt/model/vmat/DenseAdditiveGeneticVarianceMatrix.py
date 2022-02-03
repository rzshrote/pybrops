from pybropt.model.vmat.DenseGeneticVarianceMatrix import DenseGeneticVarianceMatrix
from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix

# TODO: implement me
class DenseAdditiveGeneticVarianceMatrix(DenseGeneticVarianceMatrix,AdditiveGeneticVarianceMatrix):
    """docstring for DenseAdditiveGeneticVarianceMatrix."""

    def __init__(self, arg):
        super(DenseAdditiveGeneticVarianceMatrix, self).__init__()
        self.arg = arg
