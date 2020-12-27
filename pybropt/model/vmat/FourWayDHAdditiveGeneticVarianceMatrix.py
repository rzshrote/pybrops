from . import AdditiveGeneticVarianceMatrix

class FourWayDHAdditiveGeneticVarianceMatrix(AdditiveGeneticVarianceMatrix):
    """docstring for FourWayDHAdditiveGeneticVarianceMatrix."""

    def __init__(self, arg):
        super(FourWayDHAdditiveGeneticVarianceMatrix, self).__init__()
        self.arg = arg
