from . import DenseEstimatedBreedingValueMatrix

class DensePhenotypicEstimatedBreedingValueMatrix(DenseEstimatedBreedingValueMatrix):
    """docstring for PhenotypicEstimatedBreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs: dict):
        super(DensePhenotypicEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################



################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhenotypicEstimatedBreedingValueMatrix(v):
    return isinstance(v, DensePhenotypicEstimatedBreedingValueMatrix)

def check_is_DensePhenotypicEstimatedBreedingValueMatrix(v, varname):
    if not isinstance(v, DensePhenotypicEstimatedBreedingValueMatrix):
        raise TypeError("'%s' must be a DensePhenotypicEstimatedBreedingValueMatrix." % varname)

def cond_check_is_DensePhenotypicEstimatedBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DensePhenotypicEstimatedBreedingValueMatrix(v, varname)
