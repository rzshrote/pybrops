import numpy

from pybropt.core.mat import DenseTaxaTraitMatrix
from . import BreedingValueMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim

class DenseBreedingValueMatrix(DenseTaxaTraitMatrix,BreedingValueMatrix):
    """Dense breeding value matrix implementation."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, trait = None, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        mat : numpy.ndarray
            A float64 matrix of breeding values of shape (n, t).
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.

        """
        super(DenseBreedingValueMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = trait,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################# Breeding Value Data ##################
    def mat():
        doc = "Raw matrix property."
        def fget(self):
            """Get raw matrix"""
            return self._mat
        def fset(self, value):
            """Set raw matrix"""
            check_is_ndarray(value, "mat")
            check_ndarray_ndim(value, "mat", 2)
            self._mat = value
        def fdel(self):
            """Delete raw matrix"""
            del self._mat
        return locals()
    mat = property(**mat())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############## Matrix summary statistics ###############
    def targmax(self):
        """
        Return indices of the maximum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing indices of maximum values
            along the trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.argmax(axis = 0)    # get argument maximum
        return out

    def targmin(self):
        """
        Return indices of the minimum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing indices of minimum values
            along the trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.argmin(axis = 0)    # get argument minimum
        return out

    def tmax(self):
        """
        Return the maximum along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.max(axis = 0)   # get maximum
        return out

    def tmean(self):
        """
        Return the mean along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.mean(axis = 0)  # get mean
        return out

    def tmin(self):
        """
        Return the minimum along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing minimum values along the
            trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.min(axis = 0)   # get minimum
        return out

    def trange(self):
        """
        Return the range along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = numpy.ptp(self._mat, axis = 0)    # get range
        return out

    def tstd(self):
        """
        Return the standard deviation along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing standard deviation values along
            the trait axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.std(axis = 0)   # get standard deviation
        return out

    def tvar(self):
        """
        Return the variance along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        out = self._mat.var(axis = 0)   # get variance
        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DenseBreedingValueMatrix(v):
    return isinstance(v, DenseBreedingValueMatrix)

def check_is_DenseBreedingValueMatrix(v, varname):
    if not isinstance(v, DenseBreedingValueMatrix):
        raise TypeError("'%s' must be a DenseBreedingValueMatrix." % varname)

def cond_check_is_DenseBreedingValueMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DenseBreedingValueMatrix(v, varname)
