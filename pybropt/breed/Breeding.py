class Breeding:
    """docstring for Breeding."""

    def __init__(self):
        """
        Constructor for Breeding class. Does absolutely nothing.
        """
        pass

    def objfn(self, sel, *args, **kwargs):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            Vector of selection data.
        *args : tuple
            Additional arguments for this function.
        **kwargs : dict
            Additional arguments for this function.

        Returns
        -------
        error : NotImplementedError
            This function raises an error on execution.
        """
        raise NotImplementedError("This method is not implemented.")

    def objfn_vec(self, sel, *args, **kwargs):
        """
        Breeding method objective function. Implement this in derived classes.

        Parameters
        ----------
        sel : numpy.ndarray
            A 2D matrix of selection data.
        *args : tuple
            Additional arguments for this function.
        **kwargs : dict
            Additional arguments for this function.

        Returns
        -------
        error : NotImplementedError
            This function raises an error on execution.
        """
        raise NotImplementedError("This method is not implemented.")

    @classmethod
    def optimize(self, algorithm):
        """
        """
        raise NotImplementedError("This method is not implemented.")

    @classmethod
    def simulate(self, algorithm):
        """
        """
        raise NotImplementedError("This method is not implemented.")
