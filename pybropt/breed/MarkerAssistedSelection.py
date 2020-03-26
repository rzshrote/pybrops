import breed

class MarkerAssistedSelection(breed.MolecularBreeding):
    """docstring for MarkerAssistedSelection."""

    def __init__(self, population, cross, method = "MAS"):
        super(MarkerAssistedSelection, self).__init__(population, cross, method)
