import numpy

class SetChromosomeGenerator(ChromosomeGenerator):
    """docstring for SetChromosomeGenerator."""

    def __init__(self, k, setspace, rng, **kwargs):
        super(SetChromosomeGenerator, self).__init__(**kwargs)
        self.k = k
        self.setspace = setspace
        self.rng = rng

    def generate(self, k = None, setspace = None, **kwargs):
        if k is None:
            k = self.k
        if setspace is None:
            setspace = self.setspace

        # sample without replacement
        chrom = self.rng.choice(setspace, k, False)

        return chrom
