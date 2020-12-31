from . import Matrix

class SortableMatrix(Matrix):
    """docstring for SortableMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **vargs):
        """
        SortableMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(SortableMatrix, self).__init__(**vargs)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################### Sorting Methods ####################
    def lexsort(self, keys, axis):
        raise NotImplementedError("method is abstract")

    def reorder(self, indices, axis):
        raise NotImplementedError("method is abstract")

    def sort(self, keys, axis):
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_SortableMatrix(v):
    return isinstance(v, SortableMatrix)

def check_is_SortableMatrix(v, varname):
    if not isinstance(v, SortableMatrix):
        raise TypeError("'%s' must be a SortableMatrix." % varname)

def cond_check_is_SortableMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SortableMatrix(v, varname)
