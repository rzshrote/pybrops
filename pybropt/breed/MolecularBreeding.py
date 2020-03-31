from . import Breeding

class MolecularBreeding(Breeding):
    """docstring for MolecularBreeding."""

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "Molecular"):
        super(MolecularBreeding, self).__init__(population, cross, method)
