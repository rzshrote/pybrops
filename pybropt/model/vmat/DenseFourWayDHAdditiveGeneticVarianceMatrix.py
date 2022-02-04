from pybropt.model.vmat.DenseAdditiveGeneticVarianceMatrix import DenseAdditiveGeneticVarianceMatrix

class DenseFourWayDHAdditiveGeneticVarianceMatrix(DenseAdditiveGeneticVarianceMatrix):
    """docstring for DenseFourWayDHAdditiveGeneticVarianceMatrix."""

    def __init__(self, arg):
        super(DenseFourWayDHAdditiveGeneticVarianceMatrix, self).__init__()
        self.arg = arg
