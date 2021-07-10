class BreedingValueProtocol:
    """docstring for BreedingValueProtocol."""

    def __init__(self, **kwargs):
        super(BreedingValueProtocol, self).__init__()

    def estimate(self, pt, gmat, **kwargs):
        raise NotImplementedError("method is abstract")
