from . import DenseBreedingValueMatrix

# TODO: add standard errors for this class; this could be used for two-stage estimation
class DenseGenomicEstimatedBreedingValueMatrix(DenseBreedingValueMatrix):
    """docstring for DenseGenomicEstimatedBreedingValueMatrix."""

    def __init__(self, mat, taxa = None, taxa_grp = None, trait = None, **kwargs):
        """
        Constructor for the concrete class DenseGenomicEstimatedBreedingValueMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
        taxa : numpy.ndarray
        taxa_grp : numpy.ndarray
        trait : numpy.ndarray
        **kwargs : dict
            Additional keyword arguments.
        """
        super(DenseGenomicEstimatedBreedingValueMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseGenomicEstimatedBreedingValueMatrix(v):
    """Return whether an object is a DenseGenomicEstimatedBreedingValueMatrix or not"""
    return isinstance(v, DenseGenomicEstimatedBreedingValueMatrix)

def check_is_DenseGenomicEstimatedBreedingValueMatrix(v, varname):
    """Raise TypeError if object is not a DenseGenomicEstimatedBreedingValueMatrix"""
    if not isinstance(v, DenseGenomicEstimatedBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseGenomicEstimatedBreedingValueMatrix." % varname)

def cond_check_is_DenseGenomicEstimatedBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a DenseGenomicEstimatedBreedingValueMatrix"""
    if cond(v):
        check_is_DenseGenomicEstimatedBreedingValueMatrix(v, varname)
