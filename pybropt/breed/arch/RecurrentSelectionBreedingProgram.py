from . import BreedingProgram

class RecurrentSelectionBreedingProgram(BreedingProgram):
    """docstring for RecurrentSelectionBreedingProgram."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(RecurrentSelectionBreedingProgram, self).__init__(
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############# Generation number properties #############
    def t_cur():
        doc = "Current generation number of the BreedingNode."
        def fget(self):
            return self._t_cur
        def fset(self, value):
            self._t_cur = value
        def fdel(self):
            del self._t_cur
        return locals()
    t_cur = property(**t_cur())

    def t_max():
        doc = "Maximum generation number of the BreedingNode."
        def fget(self):
            return self._t_max
        def fset(self, value):
            self._t_max = value
        def fdel(self):
            del self._t_max
        return locals()
    t_max = property(**t_max())

    ################ Population properties #################
    def geno():
        doc = "Main breeding population of the BreedingNode."
        def fget(self):
            return self._pop
        def fset(self, value):
            self._pop = value
        def fdel(self):
            del self._pop
        return locals()
    geno = property(**geno())

    ############## Breeding value properties ###############
    def bval():
        doc = "Estimated breeding values for the main breeding population of the BreedingNode."
        def fget(self):
            return self._bval
        def fset(self, value):
            self._bval = value
        def fdel(self):
            del self._bval
        return locals()
    bval = property(**bval())

    ############### Genomic model properties ###############
    def gmod():
        doc = "Estimated genomic model for the main breeding population of the BreedingNode."
        def fget(self):
            return self._gmod
        def fset(self, value):
            self._gmod = value
        def fdel(self):
            del self._gmod
        return locals()
    gmod = property(**gmod())

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
            return self._mateop
        def fset(self, value):
            self._mateop = value
        def fdel(self):
            del self._mateop
        return locals()
    mateop = property(**mateop())

    def evalop():
        doc = "Evaluation operator."
        def fget(self):
            return self._evalop
        def fset(self, value):
            self._evalop = value
        def fdel(self):
            del self._evalop
        return locals()
    evalop = property(**evalop())

    def intgop():
        doc = "Integration operator."
        def fget(self):
            return self._intgop
        def fset(self, value):
            self._intgop = value
        def fdel(self):
            del self._intgop
        return locals()
    intgop = property(**intgop())

    def calop():
        doc = "Genomic model calibration operator."
        def fget(self):
            return self._calop
        def fset(self, value):
            self._calop = value
        def fdel(self):
            del self._calop
        return locals()
    calop = property(**calop())

    def sselop():
        doc = "Survivor selection operator."
        def fget(self):
            return self._sselop
        def fset(self, value):
            self._sselop = value
        def fdel(self):
            del self._sselop
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

def check_is_RecurrentSelectionBreedingProgram(v, vname):
    if not isinstance(v, RecurrentSelectionBreedingProgram):
        raise TypeError("variable '{0}' must be a RecurrentSelectionBreedingProgram".format(vname))

def cond_check_is_BreedingProgram(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_RecurrentSelectionBreedingProgram(v, vname)
