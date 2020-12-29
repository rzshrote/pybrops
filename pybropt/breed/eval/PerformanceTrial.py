class PerformanceTrial:
    """docstring for PerformanceTrial."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PerformanceTrial, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def evaluate(self, **kwargs):
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PerformanceTrial(v):
    return isinstance(v, PerformanceTrial)

def check_is_PerformanceTrial(v, varname):
    if not isinstance(v, PerformanceTrial):
        raise TypeError("'%s' must be a PerformanceTrial." % varname)

def cond_check_is_PerformanceTrial(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PerformanceTrial(v, varname)
