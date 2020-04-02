from . import GenomicMating

class SPstdA(GenomicMating):
    """
    Sum of Progeny Standard Deviations of Additive effects (SPstdA) objective
    function class.
    """

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "SPstdA"):
        super(SPstdA, self).__init__(population, cross, method)

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################
    def optimize(self, algorithm, objcoeff = None, minimizing = True, **kwargs):
        sel = super(SPstdA, self).optimize(
            algorithm = algorithm,
            objcoeff = objcoeff,
            minimizing = minimizing,
            **kwargs
        )
        return sel

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def objfn_mat(sel, geno, coeff, varAfn):
        raise NotImplementedError
