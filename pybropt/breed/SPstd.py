class SPstd(GenomicMating):
    """
    Sum of Progeny STandard Deviations (SPstd) objective function class.
    """

    def __init__(self, population, cross):
        super(SPstd, self).__init__(population, cross)
