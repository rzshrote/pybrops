from . import Breeding

class PhenotypicBreeding(Breeding):
    """docstring for PhenotypicBreeding."""

    ############################################################################
    ######################### Reserved Object Methods ##########################
    ############################################################################
    def __init__(self, population, cross, method = "Phenotypic"):
        super(PhenotypicBreeding, self).__init__(population, cross, method)