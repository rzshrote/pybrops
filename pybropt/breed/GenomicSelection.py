import breed

class GenomicSelection(breed.MolecularBreeding):
    """docstring for GenomicSelection."""

    def __init__(self, population, cross, method = "GS"):
        super(GenomicSelection, self).__init__(population, cross, method)
