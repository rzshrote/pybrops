class GeneticAlgorithm:
    """docstring for GeneticAlgorithm."""

    def __init__(self, **kwargs):
        super(GeneticAlgorithm, self).__init__()

    def mu():
        doc = "The mu property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    mu = property(**mu())

    def lamb():
        doc = "The lamb property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    lamb = property(**lamb())

    def popgen():
        doc = "The popgen property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    popgen = property(**popgen())

    def pselop():
        doc = "The pselop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    pselop = property(**pselop())

    def cxop():
        doc = "The cxop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    cxop = property(**cxop())

    def mutop():
        doc = "The mutop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    mutop = property(**mutop())

    def memeop():
        doc = "The memeop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    memeop = property(**memeop())

    def sselop():
        doc = "The sselop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    sselop = property(**sselop())

    def evolve(self, t, lbook, **kwargs):
        raise NotImplementedError("method is abstract")
