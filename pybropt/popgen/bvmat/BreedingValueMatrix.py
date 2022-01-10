from pybropt.core.mat import TaxaTraitMatrix
from pybropt.core.io.HDF5InputOutput import HDF5InputOutput

class BreedingValueMatrix(TaxaTraitMatrix,HDF5InputOutput):
    """
    The BreedingValueMatrix class represents a Multivariate Breeding Value.

    All elements within a BreedingValueMatrix are mean-centered and scaled to
    unit variance for each trait.
                                     X - μ
                                BV = ─────
                                       σ
        Where:
            BV is the breeding value.
            X is the phenotype value.
            μ is the mean (location) for X.
            σ is the standard deviation (scale) for X.

    Phenotype values can be reconstituted using:
        X = σBV + μ
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        BreedingValueMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(BreedingValueMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def location():
        doc = "Mean of the phenotype values used to calculate breeding values"
        def fget(self):
            """Get the mean of the phenotype values used to calculate breeding values"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the mean of the phenotype values used to calculate breeding values"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the mean of the phenotype values used to calculate breeding values"""
            raise NotImplementedError("method is abstract")
        return locals()
    location = property(**location())

    def scale():
        doc = "Standard deviation of the phenotype values used to calculate breeding values"
        def fget(self):
            """Get the standard deviation of the phenotype values used to calculate breeding values"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the standard deviation of the phenotype values used to calculate breeding values"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the standard deviation of the phenotype values used to calculate breeding values"""
            raise NotImplementedError("method is abstract")
        return locals()
    scale = property(**scale())

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

    def tmax(self, descale):
        """
        Return the maximum along the trait axis.

        Parameters
        ----------
        descale : boolean
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tmean(self, descale):
        """
        Return the mean along the trait axis.

        Parameters
        ----------
        descale : boolean
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing maximum values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tmin(self, descale):
        """
        Return the minimum along the trait axis.

        Parameters
        ----------
        descale : boolean
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape (t,) containing minimum values along the
            trait axis.
            Where:
                't' is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def trange(self, descale):
        """
        Return the range along the trait axis.

        Parameters
        ----------
        descale : boolean
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tstd(self, descale):
        """
        Return the standard deviation along the trait axis.

        Parameters
        ----------
        descale : boolean
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing standard deviation values along
            the trait axis.
            Where:
                't' is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tvar(self, descale):
        """
        Return the variance along the trait axis.

        Parameters
        ----------
        descale : boolean
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (t,) containing variance values along the trait
            axis.
            Where:
                't' is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def descale(self):
        """
        Transform values within the BreedingValueMatrix back to their de-scaled
        and de-centered values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape (n,t) containing de-scaled and de-centered values.
            Where:
                'n' is the number of taxa.
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
