class PhenotypingProtocol:
    """docstring for PhenotypingProtocol."""

    def __init__(self, **kwargs):
        super(PhenotypingProtocol, self).__init__()

    def phenotype(self, pgmat, gpmod, **kwargs):
        raise NotImplementedError("method is abstract")
