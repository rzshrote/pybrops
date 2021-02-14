class CrossoverOperator:
    """docstring for CrossoverOperator."""

    def __init__(self, **kwargs):
        super(CrossoverOperator, self).__init__()

    def nparent():
        doc = "The nparent property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    nparent = property(**nparent())

    def crossover(self, parents, **kwargs):
        """
        Parameters
        ----------
        parents : tuple or list_like

        Returns
        -------
        offspring : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")
