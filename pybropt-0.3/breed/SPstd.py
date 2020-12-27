from . import GenomicMating

class SPstd(GenomicMating):
    """
    Sum of Progeny STandard Deviations (SPstd) objective function class.
    """

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "SPstd"):
        super(SPstd, self).__init__(population, cross, method)

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def optimize(self, algorithm, objcoeff = None, minimizing = True, **kwargs):
        sel = super(SPstd, self).optimize(
            algorithm = algorithm,
            objcoeff = objcoeff,
            minimizing = minimizing,
            **kwargs
        )
        return sel
