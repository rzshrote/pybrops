from pybropt.base.mat import TaxaMatrix

class BreedingValueMatrix(TaxaMatrix):
    """docstring for BreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(BreedingValueMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Breeding Value Data ##################
    def trait():
        doc = "The trait property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    trait = property(**trait())

    def ntrait():
        doc = "The ntrait property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    ntrait = property(**ntrait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################