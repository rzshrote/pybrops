from pybropt.model.vmat.AdditiveGeneticVarianceMatrix import AdditiveGeneticVarianceMatrix

class DihybridDHAdditiveGeneticVarianceMatrix(AdditiveGeneticVarianceMatrix):
    """docstring for DihybridDHAdditiveGeneticVarianceMatrix."""

    def __init__(self, arg):
        super(DihybridDHAdditiveGeneticVarianceMatrix, self).__init__()
        self.arg = arg
