class NoChromosomeDecoder(ChromosomeDecoder):
    """docstring for NoChromosomeDecoder."""

    def __init__(self, **kwargs):
        super(NoChromosomeDecoder, self).__init__(**kwargs)

    def decode(self, chrom, **kwargs):
        """
        Parameters
        ----------
        chrom : numpy.ndarray

        Returns
        -------
        soln : list_like
        """
        return chrom
