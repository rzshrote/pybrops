from . import BreedingNode

class GermplasmBank(BreedingNode):
    """docstring for GermplasmBank."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(GermplasmBank, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_GermplasmBank(v):
    return isinstance(v, GermplasmBank)

def check_is_GermplasmBank(v, varname):
    if not isinstance(v, GermplasmBank):
        raise TypeError("'%s' must be a GermplasmBank." % varname)

def cond_check_is_GermplasmBank(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GermplasmBank(v, varname)
