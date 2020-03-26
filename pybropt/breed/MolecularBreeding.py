import breed

class MolecularBreeding(breed.Breeding):
    """docstring for MolecularBreeding."""

    ############################################################################
    ############################# Reserved methods #############################
    ############################################################################
    def __init__(self, population, cross, method = "Molecular"):
        super(MolecularBreeding, self).__init__(population, cross, method)
