class PopulationGenerator:
    """docstring for PopulationGenerator."""

    def __init__(self, **kwargs: dict):
        super(PopulationGenerator, self).__init__()

    def generate(self, nchrom, **kwargs: dict):
        raise NotImplementedError("method is abstract")
