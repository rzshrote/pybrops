from . import DenseBreedingValueMatrix

# TODO: add standard errors for this class; this could be used for two-stage estimation
class DenseEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseEstimatedBreedingValueMatrix."""

    def __init__(self, mat, taxa = None, taxa_grp = None, trait = None, **kwargs):
        super(DenseEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseEstimatedBreedingValueMatrix(v):
    return isinstance(v, DenseEstimatedBreedingValueMatrix)

def check_is_DenseEstimatedBreedingValueMatrix(v, varname):
    if not isinstance(v, DenseEstimatedBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseEstimatedBreedingValueMatrix." % varname)

def cond_check_is_DenseEstimatedBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseEstimatedBreedingValueMatrix(v, varname)
