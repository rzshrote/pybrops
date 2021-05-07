class MemeticOperator:
    """docstring for MemeticOperator."""

    def __init__(self, **kwargs):
        super(MemeticOperator, self).__init__()

    def evolve(self, objfn, pop, score, **kwargs):
        """
        Parameters
        ----------
        objfn : function
        pop : list of numpy.ndarray
        score : numpy.ndarray
        **kwargs :
            Additional arguments to pass to objfn.

        Returns
        -------
        new_pop : list of numpy.ndarray
        new_score : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")
