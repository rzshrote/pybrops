class ChromosomeDecoder:
    """docstring for ChromosomeDecoder."""

    def __init__(self, **kwargs):
        super(ChromosomeDecoder, self).__init__()

    def decode(self, chrom, **kwargs):
        """
        Parameters
        ----------
        chrom : numpy.ndarray

        Returns
        -------
        soln : list_like
        """
        raise NotImplementedError("method is abstract")
