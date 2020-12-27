from . import GeneticVarianceMatrix

class AdditiveGeneticVarianceMatrix(GeneticVarianceMatrix):
    """docstring for AdditiveGeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(AdditiveGeneticVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Variance Matrix Data #################
    def var_A():
        doc = "The var_A property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    var_A = property(**var_A())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
