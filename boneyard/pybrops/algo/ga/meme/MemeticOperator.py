class MemeticOperator:
    """docstring for MemeticOperator."""

    def __init__(self, **kwargs: dict):
        super(MemeticOperator, self).__init__()

    def evolve(self, objfn, pop, score, **kwargs: dict):
        """
        Parameters
        ----------
        objfn : function
        pop : list of numpy.ndarray
        score : numpy.ndarray
        kwargs : dict
            Additional arguments to pass to objfn.

        Returns
        -------
        new_pop : list of numpy.ndarray
        new_score : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")
