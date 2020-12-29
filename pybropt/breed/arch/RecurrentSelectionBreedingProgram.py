from . import BreedingProgram

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """docstring for RecurrentSelectionBreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, arg):
        super(RecurrentSelectionBreedingProgram, self).__init__()
        self.arg = arg

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############# Generation number properties #############
    def gen_cur():
        doc = "Current generation number of the BreedingNode."
        def fget(self):
            return self._gen_cur
        def fset(self, value):
            self._gen_cur = value
        def fdel(self):
            del self._gen_cur
        return locals()
    gen_cur = property(**gen_cur())

    def gen_max():
        doc = "Maximum generation number of the BreedingNode."
        def fget(self):
            return self._gen_max
        def fset(self, value):
            self._gen_max = value
        def fdel(self):
            del self._gen_max
        return locals()
    gen_max = property(**gen_max())

    ################ Population properties #################
    def pop():
        doc = "Main breeding population of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    pop = property(**pop())

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
    def bval():
        doc = "Estimated breeding values for the main breeding population of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    bval = property(**bval())

    def bval_true():
        doc = "True breeding values for the main breeding population of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    bval_true = property(**bval_true())

    ############### Genomic model properties ###############
    def gmod():
        doc = "Estimated genomic model for the main breeding population of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    gmod = property(**gmod())

    def gmod_true():
        doc = "True genomic model for the main breeding population of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    gmod_true = property(**gmod_true())

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
def is_RecurrentSelectionBreedingProgram(v):
    return isinstance(v, RecurrentSelectionBreedingProgram)

def check_is_RecurrentSelectionBreedingProgram(v, varname):
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("'%s' must be a RecurrentSelectionBreedingProgram." % varname)

def cond_check_is_BreedingProgram(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_RecurrentSelectionBreedingProgram(v, varname)
