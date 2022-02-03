from pybropt.core.mat.DenseSquareTaxaMatrix import DenseSquareTaxaMatrix
from pybropt.model.vmat.GeneticVarianceMatrix import GeneticVarianceMatrix

# TODO: implement me
class DenseGeneticVarianceMatrix(DenseSquareTaxaMatrix,GeneticVarianceMatrix):
    """docstring for DenseGeneticVarianceMatrix."""

    def __init__(self, arg):
        super(DenseGeneticVarianceMatrix, self).__init__()
        self.arg = arg
        
