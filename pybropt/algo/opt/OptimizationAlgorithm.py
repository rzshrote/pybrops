class OptimizationAlgorithm:
    """docstring for OptimizationAlgorithm."""

    def __init__(self, **kwargs):
        """
        Constructor for the abstract class OptimizationAlgorithm.
        """
        super(OptimizationAlgorithm, self).__init__()

    def optimize(self, objfn, k, sspace, objwt, **kwargs):
        """
        Optimize an objective function.

        Parameters
        ----------
        objfn : callable
            Objective function which to optimize.
        k : int
            Number of decision variables in the search space.
            A vector is formed as sspace^k
        sspace : numpy.ndarray
            Search space that the OptimizationAlgorithm searches in.
        objwt : numpy.ndarray
            Weight(s) applied to output(s) from the objfn.

        Returns
        -------
        out : tuple
            A tuple of length 3 (soln, decn, misc)
        """
        raise NotImplementedError("method is abstract")
