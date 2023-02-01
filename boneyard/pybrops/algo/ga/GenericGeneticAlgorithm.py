from pybrops.algo.opt.ga.GeneticAlgorithm import GeneticAlgorithm

class GenericGeneticAlgorithm(GeneticAlgorithm):
    """
    docstring for GenericGeneticAlgorithm.
    """
    def __init__(self, mu, lamb, ngen, **kwargs: dict):
        """
        Constructor for GenericGeneticAlgorithm.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(GenericGeneticAlgorithm, self).__init__(**kwargs)
        self.mu = mu
        self.lamb = lamb

    def mu():
        doc = "The mu property."
        def fget(self):
            """Get value for mu."""
            return self._mu
        def fset(self, value):
            """Set value for mu."""
            self._mu = value
        def fdel(self):
            """Delete value for mu."""
            del self._mu
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    mu = property(**mu())

    def lamb():
        doc = "The lamb property."
        def fget(self):
            """Get value for lamb."""
            return self._lamb
        def fset(self, value):
            """Set value for lamb."""
            self._lamb = value
        def fdel(self):
            """Delete value for lamb."""
            del self._lamb
        return {"fget":fget, "fset":fset, "fdel":fdel, "doc":doc}
    lamb = property(**lamb())

    def evolve(self, t):
        """
        Evolve population for t generations.
        
        Parameters
        ----------
        t : int
            Number of generations for which to evolve.
        
        Returns
        -------
        out : outtype
            outdesc
        """
        pass