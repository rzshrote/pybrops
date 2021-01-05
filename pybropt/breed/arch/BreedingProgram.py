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

    ############## Breeding value properties ###############
    def bval_queue():
        doc = "The bval_queue property."
        def fget(self):
            return self._bval_queue
        def fset(self, value):
            self._bval_queue = value
        def fdel(self):
            del self._bval_queue
        return locals()
    bval_queue = property(**bval_queue())

    def bval_kw():
        doc = "The bval_kw property."
        def fget(self):
            return self._bval_kw
        def fset(self, value):
            self._bval_kw = value
        def fdel(self):
            del self._bval_kw
        return locals()
    bval_kw = property(**bval_kw())

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
    def pselect(gen_cur, gen_max, pop, bval, bval_true, gmod, gmod_true,):
        pass


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
