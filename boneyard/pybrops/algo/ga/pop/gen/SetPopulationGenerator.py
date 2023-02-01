from . import PopulationGenerator

from pybrob.algo.ga.chrom.gen import SetChromosomeGenerator

class SetPopulationGenerator(PopulationGenerator):
    """docstring for SetPopulationGenerator."""

    def __init__(self, n, k, setspace, rng, **kwargs: dict):
        super(SetPopulationGenerator, self).__init__(**kwargs)
        self.n = n
        self.k = k
        self.setspace = setspace
        self.rng = rng
        self.chrgenop = SetChromosomeGenerator(
            k = self.k,
            setspace = self.setspace,
            rng = self.rng
        )

    def generate(self, n = None, k = None, setspace = None, **kwargs: dict):
        if n is None:
            n = self.n
        if k is None:
            k = self.k
        if setspace is None:
            setspace = self.setspace

        pop = []
        for i in range(n):
            pop.append(
                self.chrgenop.generate(
                    k = k,
                    setspace = setspace,
                    **kwargs
                )
            )

        return pop
