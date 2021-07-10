class GenotypingProtocol:
    """docstring for GenotypingProtocol."""

    def __init__(self, **kwargs):
        super(GenotypingProtocol, self).__init__()

    def genotype(self, pgmat, **kwargs):
        raise NotImplementedError("method is abstract")
