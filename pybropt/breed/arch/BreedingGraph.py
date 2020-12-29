

class BreedingGraph:
    """docstring for BreedingGraph."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BreedingGraph, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Graph properties ###################
    def graph():
        doc = "Access to the underlying graph structure."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    graph = property(**graph())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingGraph(v):
    return isinstance(v, BreedingGraph)

def check_is_BreedingGraph(v, varname):
    if not isinstance(v, BreedingGraph):
        raise TypeError("'%s' must be a BreedingGraph." % varname)

def cond_check_is_BreedingGraph(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingGraph(v, varname)
