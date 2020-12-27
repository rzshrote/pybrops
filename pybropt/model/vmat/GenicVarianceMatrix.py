from . import VarianceMatrix

class GenicVarianceMatrix(VarianceMatrix):
    """docstring for GenicVarianceMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GenicVarianceMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Variance Matrix Data #################
    def var_a():
        doc = "The var_a property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    var_a = property(**var_a())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
