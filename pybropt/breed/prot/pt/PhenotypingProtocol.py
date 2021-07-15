class PhenotypingProtocol:
    """docstring for PhenotypingProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(PhenotypingProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def phenotype(self, pgmat, gpmod, **kwargs):
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_PhenotypingProtocol(v):
    return isinstance(v, PhenotypingProtocol)

def check_is_PhenotypingProtocol(v, vname):
    if not isinstance(v, PhenotypingProtocol):
        raise TypeError("variable '{0}' must be a PhenotypingProtocol".format(vname))

def cond_check_is_PhenotypingProtocol(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_PhenotypingProtocol(v, vname)
