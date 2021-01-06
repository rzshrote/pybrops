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
    def t_cur():
        doc = "Current generation number of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    t_cur = property(**t_cur())

    def t_max():
        doc = "Maximum generation number of the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    t_max = property(**t_max())

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
        doc = "Dictionary of breeding values for the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    bval = property(**bval())

    ############### Genomic model properties ###############
    def gmod():
        doc = "Dictionary of genomic models for the BreedingNode."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    gmod = property(**gmod())

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
