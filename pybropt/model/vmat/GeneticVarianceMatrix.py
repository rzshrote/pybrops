from . import VarianceMatrix

class GeneticVarianceMatrix(VarianceMatrix):
    """docstring for GeneticVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GeneticVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Variance Matrix Data #################
    def vmat_G():
        doc = "The vmat_G property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vmat_G = property(**vmat_G())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
