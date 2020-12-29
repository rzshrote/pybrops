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

    ################ Population properties #################
    def pop_queue():
        doc = "Queue breeding population of the BreedingProgram."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    pop_queue = property(**pop_queue())

    def pop_kw():
        doc = "Keyword breeding population of the BreedingProgram."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    pop_kw = property(**pop_kw())

    ######### Breeding program operator properties #########
    def pselop():
        doc = "The pselop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    pselop = property(**pselop())

    def mateop():
        doc = "The mateop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    mateop = property(**mateop())

    def evalop():
        doc = "The evalop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    evalop = property(**evalop())

    def calop():
        doc = "The calop property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    calop = property(**calop())

    def sselop():
        doc = "The sselop property."
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingProgram(v):
    return isinstance(v, BreedingProgram)

def check_is_BreedingProgram(v, varname):
    if not isinstance(v, BreedingProgram):
        raise TypeError("'%s' must be a BreedingProgram." % varname)

def cond_check_is_BreedingProgram(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingProgram(v, varname)
