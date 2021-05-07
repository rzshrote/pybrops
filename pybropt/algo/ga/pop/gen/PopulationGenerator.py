class PopulationGenerator:
    """docstring for PopulationGenerator."""

    def __init__(self, **kwargs):
        super(PopulationGenerator, self).__init__()

    def generate(self, nchrom, **kwargs):
        raise NotImplementedError("method is abstract")
