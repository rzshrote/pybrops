from . import MolecularBreeding

class MarkerAssistedSelection(MolecularBreeding):
    """docstring for MarkerAssistedSelection."""

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "MAS"):
        super(MarkerAssistedSelection, self).__init__(population, cross, method)
