from . import MolecularBreeding

class GenomicSelection(MolecularBreeding):
    """docstring for GenomicSelection."""

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "GS"):
        super(GenomicSelection, self).__init__(population, cross, method)