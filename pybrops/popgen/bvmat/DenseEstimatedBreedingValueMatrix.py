from pybrops.popgen.bvmat.DenseBreedingValueMatrix import DenseBreedingValueMatrix

# TODO: add standard errors for this class; this could be used for two-stage estimation
class DenseEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseEstimatedBreedingValueMatrix."""

    def __init__(self, mat, location, scale, taxa = None, taxa_grp = None, trait = None, **kwargs):
        """
        Constructor for the concrete class DenseEstimatedBreedingValueMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape (n, t).
        location : numpy.ndarray
            An array of breeding value locations of shape (t,).
        scale : numpy.ndarray
            An array of breeding value scales of shape (t,).
        taxa : numpy.ndarray
            An array of taxa names.
        taxa_grp : numpy.ndarray
            An array of taxa groups.
        trait : numpy.ndarray
            An array of trait names.
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DenseEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            location = location,
            scale = scale,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseEstimatedBreedingValueMatrix(v):
    """Return whether an object is a DenseEstimatedBreedingValueMatrix or not"""
    return isinstance(v, DenseEstimatedBreedingValueMatrix)

def check_is_DenseEstimatedBreedingValueMatrix(v, varname):
    """Raise TypeError if object is not a DenseEstimatedBreedingValueMatrix"""
    if not isinstance(v, DenseEstimatedBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseEstimatedBreedingValueMatrix." % varname)

def cond_check_is_DenseEstimatedBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a DenseEstimatedBreedingValueMatrix"""
    if cond(v):
        check_is_DenseEstimatedBreedingValueMatrix(v, varname)
