from . import BreedingNode

class BreedingProgram(BreedingNode):
    """docstring for BreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BreedingProgram, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Selection age properties ###############
    def a_max():
        doc = "The a_max property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    a_max = property(**a_max())

    ######### Breeding program operator properties #########
    def pselop():
        doc = "Parental selection operator."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    pselop = property(**pselop())

    def mateop():
        doc = "Mating operator."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    mateop = property(**mateop())

    def evalop():
        doc = "Evaluation operator."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    evalop = property(**evalop())

    def intgop():
        doc = "Integration operator."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    intgop = property(**intgop())

    def calop():
        doc = "Genomic model calibration operator."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    calop = property(**calop())

    def sselop():
        doc = "Survivor selection operator."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    sselop = property(**sselop())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################ Whole breeding program ################
    def evolve(self, **kwargs):
        raise NotImplementedError("method is abstract")

################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingProgram(v):
    return isinstance(v, BreedingProgram)

def check_is_BreedingProgram(v, vname):
    if not isinstance(v, BreedingProgram):
        raise TypeError("variable '{0}' must be a BreedingProgram".format(vname))

def cond_check_is_BreedingProgram(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingProgram(v, vname)
