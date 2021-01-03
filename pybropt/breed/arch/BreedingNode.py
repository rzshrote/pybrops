

class BreedingNode:
    """docstring for BreedingNode."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BreedingNode, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############# Generation number properties #############
    def gen_cur():
        doc = "Current generation number of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    gen_cur = property(**gen_cur())

    def gen_max():
        doc = "Maximum generation number of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
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

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingNode(v):
    return isinstance(v, BreedingNode)

def check_is_BreedingNode(v, vname):
    if not isinstance(v, BreedingNode):
        raise TypeError("variable '{0}' must be a BreedingNode".format(vname))

def cond_check_is_BreedingNode(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingNode(v, vname)
