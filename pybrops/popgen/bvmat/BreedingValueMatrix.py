"""
Module defining basal matrix interfaces and associated error checking routines
for breeding value matrices.
"""

from typing import Any

import numpy
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat.TaxaTraitMatrix import TaxaTraitMatrix

class BreedingValueMatrix(TaxaTraitMatrix,HDF5InputOutput):
    """
    The BreedingValueMatrix class represents a Multivariate Breeding Value.

    Notes
    -----
    All elements within a BreedingValueMatrix are mean-centered and scaled to
    unit variance for each trait.

    .. math::
        BV = \\frac{X - \\mu}{\\sigma}

    Where:

    - :math:`BV` is the breeding value.
    - :math:`X` is the phenotype value.
    - :math:`\\mu` is the mean (location) for :math:`X`.
    - :math:`\\sigma` is the standard deviation (scale) for :math:`X`.

    Phenotype values can be reconstituted using:

    .. math::
        X = \\sigma BV + \\mu
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the abstract class BreedingValueMatrix.

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
    @property
    def location(self) -> Any:
        """Mean of the phenotype values used to calculate breeding values."""
        raise NotImplementedError("property is abstract")
    @location.setter
    def location(self, value: Any) -> None:
        """Set the mean of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")
    @location.deleter
    def location(self) -> None:
        """Delete the mean of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")
    
    @property
    def scale(self) -> Any:
        """Standard deviation of the phenotype values used to calculate breeding values."""
        raise NotImplementedError("property is abstract")
    @scale.setter
    def scale(self, value: Any) -> None:
        """Set the standard deviation of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")
    @scale.deleter
    def scale(self) -> None:
        """Delete the standard deviation of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############## Matrix summary statistics ###############
    def targmax(self) -> numpy.ndarray:
        """
        Return indices of the maximum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of maximum
            values along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def targmin(self) -> numpy.ndarray:
        """
        Return indices of the minimum values along the trait axis.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing indices of minimum
            values along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tmax(self, descale: bool) -> numpy.ndarray:
        """
        Return the maximum along the trait axis.

        Parameters
        ----------
        descale : bool
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the trait
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tmean(self, descale: bool) -> numpy.ndarray:
        """
        Return the mean along the trait axis.

        Parameters
        ----------
        descale : bool
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing maximum values along the trait
            axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tmin(self, descale: bool) -> numpy.ndarray:
        """
        Return the minimum along the trait axis.

        Parameters
        ----------
        descale : bool
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An index array of shape ``(t,)`` containing minimum values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def trange(self, descale: bool) -> numpy.ndarray:
        """
        Return the range along the trait axis.

        Parameters
        ----------
        descale : bool
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tstd(self, descale: bool) -> numpy.ndarray:
        """
        Return the standard deviation along the trait axis.

        Parameters
        ----------
        descale : bool
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing standard deviation values
            along the trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def tvar(self, descale: bool) -> numpy.ndarray:
        """
        Return the variance along the trait axis.

        Parameters
        ----------
        descale : bool
            whether to transform results to their de-scaled values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing variance values along the
            trait axis.

            Where:

            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def descale(self) -> numpy.ndarray:
        """
        Transform values within the BreedingValueMatrix back to their de-scaled
        and de-centered values.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(n,t)`` containing de-scaled and de-centered
            values.

            Where:

            - ``n`` is the number of taxa.
            - ``t`` is the number of traits.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingValueMatrix(v: object) -> bool:
    """
    Determine whether an object is a BreedingValueMatrix.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a BreedingValueMatrix object instance.
    """
    return isinstance(v, BreedingValueMatrix)

def check_is_BreedingValueMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type BreedingValueMatrix. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, BreedingValueMatrix):
        raise TypeError("variable '{0}' must be a BreedingValueMatrix".format(vname))
