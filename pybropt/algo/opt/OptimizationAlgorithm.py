class OptimizationAlgorithm:
    """docstring for OptimizationAlgorithm."""

    def __init__(self, **kwargs):
        super(OptimizationAlgorithm, self).__init__()

    def optimize(self, objfn, **kwargs):
        raise NotImplementedError("method is abstract")
