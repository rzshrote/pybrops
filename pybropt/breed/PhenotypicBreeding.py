import breed

class PhenotypicBreeding(breed.Breeding):
    """docstring for PhenotypicBreeding."""

    def __init__(self, population, cross, method = "Phenotypic"):
        super(PhenotypicBreeding, self).__init__(population, cross, method)
