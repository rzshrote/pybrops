from pybropt.core.mat import TaxaMatrix

class CoancestryMatrix(TaxaMatrix):
    """docstring for CoancestryMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(CoancestryMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genotype Data Properites ###############
    # gmat should be implemented in a sparse version of this matrix.

    ############## Coancestry Data Properites ##############
    # access using mat (inherited from Matrix)

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################## Coancestry Methods ##################
    def coancestry(*args):
        """
        Retrieve the coancestry between individuals in *args.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_CoancestryMatrix(v):
    return isinstance(v, CoancestryMatrix)

def check_is_CoancestryMatrix(v, varname):
    if not isinstance(v, CoancestryMatrix):
        raise TypeError("'%s' must be a CoancestryMatrix." % varname)

def cond_check_is_CoancestryMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_CoancestryMatrix(v, varname)
