from pybropt.core.mat import TaxaTraitMatrix
from pybropt.core.io import HDF5InputOutput

class BreedingValueMatrix(TaxaTraitMatrix,HDF5InputOutput):
    """docstring for BreedingValueMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(BreedingValueMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingValueMatrix(v):
    """Return whether an object is a BreedingValueMatrix or not"""
    return isinstance(v, BreedingValueMatrix)

def check_is_BreedingValueMatrix(v, vname):
    """Raise TypeError if object is not a BreedingValueMatrix"""
    if not isinstance(v, BreedingValueMatrix):
        raise TypeError("variable '{0}' must be a BreedingValueMatrix".format(vname))

def cond_check_is_BreedingValueMatrix(v, vname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a BreedingValueMatrix"""
    if cond(v):
        check_is_BreedingValueMatrix(v, vname)
