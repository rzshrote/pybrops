class SPstdA(GenomicMating):
    """
    Sum of Progeny Standard Deviations of Additive effects (SPstdA) objective
    function class.
    """

    def __init__(self, population, cross):
        super(SPstdA, self).__init__(population, cross)

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def objfn_mat(sel, geno, coeff, varAfn):
        raise NotImplementedError
