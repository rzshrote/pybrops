from pybropt.core.mat import Matrix

class HaplotypeMatrix(Matrix):
    """docstring for HaplotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        HaplotypeMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(HaplotypeMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def haplo():
        doc = "The haplo property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    haplo = property(**haplo())

    def ploidy():
        doc = "The ploidy property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    ploidy = property(**ploidy())

    def ntaxa():
        doc = "The ntaxa property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    ntaxa = property(**ntaxa())

    def nloci():
        doc = "The nloci property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    nloci = property(**nloci())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
