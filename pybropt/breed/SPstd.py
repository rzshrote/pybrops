from . import GenomicMating

class SPstd(GenomicMating):
    """
    Sum of Progeny STandard Deviations (SPstd) objective function class.
    """

    def __init__(self, population, cross, method = "SPstd"):
        super(SPstd, self).__init__(population, cross, method)
