from . import AdditiveGeneticVarianceMatrix

class ThreeWayDHAdditiveGeneticVarianceMatrix(AdditiveGeneticVarianceMatrix):
    """docstring for ThreeWayDHAdditiveGeneticVarianceMatrix."""

    def __init__(self, arg):
        super(ThreeWayDHAdditiveGeneticVarianceMatrix, self).__init__()
        self.arg = arg
        
