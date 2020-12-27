from pybropt.popgen.mat import TaxaMatrix

class VarianceMatrix(TaxaMatrix):
    """docstring for VarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(VarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Variance Matrix Data #################
    def vmat():
        doc = "The vmat property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    vmat = property(**vmat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
