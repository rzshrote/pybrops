class SortableMatrixInterface:
    """docstring for SortableMatrixInterface."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(SortableMatrixInterface, self).__init__(**kwargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys, axis, **kwargs):
        raise NotImplementedError("method is abstract")

    def reorder(self, indices, axis, **kwargs):
        raise NotImplementedError("method is abstract")

    def sort(self, keys, axis, **kwargs):
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_SortableMatrixInterface(v):
    return isinstance(v, SortableMatrixInterface)

def check_is_SortableMatrixInterface(v, varname):
    if not isinstance(v, SortableMatrixInterface):
        raise TypeError("'%s' must be a SortableMatrixInterface." % varname)

def cond_check_is_SortableMatrixInterface(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SortableMatrixInterface(v, varname)
