from . import MolecularBreeding

class GenomicMating(MolecularBreeding):
    """docstring for GenomicMating."""

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "GM"):
        super(GenomicMating, self).__init__(population, cross, method)
