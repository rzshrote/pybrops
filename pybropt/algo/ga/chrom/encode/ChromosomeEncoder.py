class ChromosomeEncoder:
    """docstring for ChromosomeEncoder."""

    def __init__(self, **kwargs):
        super(ChromosomeEncoder, self).__init__()

    def encode(self, soln, **kwargs):
        """
        Parameters
        ----------
        soln : list_like

        Returns
        -------
        chrom : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")
