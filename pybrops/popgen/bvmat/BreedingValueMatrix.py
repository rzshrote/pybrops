"""
Module defining basal matrix interfaces and associated error checking routines
for breeding value matrices.
"""

__all__ = [
    "BreedingValueMatrix",
    "check_is_BreedingValueMatrix",
]

from abc import ABCMeta, abstractmethod
from typing import Optional
import numpy
import pandas
from pybrops.core.io.CSVInputOutput import CSVInputOutput
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.io.PandasInputOutput import PandasInputOutput
from pybrops.core.mat.TaxaTraitMatrix import TaxaTraitMatrix

class BreedingValueMatrix(TaxaTraitMatrix,PandasInputOutput,CSVInputOutput,HDF5InputOutput,metaclass=ABCMeta):
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

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################
    @property
    @abstractmethod
    def location(self) -> object:
        """Mean of the phenotype values used to calculate breeding values."""
        raise NotImplementedError("property is abstract")
    @location.setter
    @abstractmethod
    def location(self, value: object) -> None:
        """Set the mean of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def scale(self) -> object:
        """Standard deviation of the phenotype values used to calculate breeding values."""
        raise NotImplementedError("property is abstract")
    @scale.setter
    @abstractmethod
    def scale(self, value: object) -> None:
        """Set the standard deviation of the phenotype values used to calculate breeding values"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ############## Matrix summary statistics ###############
    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    @abstractmethod
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

    #################### Export Methods ####################
    
    # inherit to_pandas (including docstrings)
    # inherit to_csv (including docstrings)
    # inherit to_hdf5 (including docstrings)

    ############################## Class Methods ###############################

    #################### Import Methods ####################
    @classmethod
    @abstractmethod
    def from_pandas(
            cls, 
            df: pandas.DataFrame, 
            **kwargs: dict
        ) -> 'BreedingValueMatrix':
        """
        Read a BreedingValueMatrix from a pandas.DataFrame.

        Parameters
        ----------
        df : pandas.DataFrame
            Pandas dataframe from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            pandas.DataFrame.

        Returns
        -------
        out : BreedingValueMatrix
            A BreedingValueMatrix read from a pandas.DataFrame.
        """
        raise NotImplementedError("class method is abstract")
    
    @classmethod
    @abstractmethod
    def from_csv(
            cls, 
            filename: str, 
            **kwargs: dict
        ) -> 'BreedingValueMatrix':
        """
        Read a BreedingValueMatrix from a CSV file.

        Parameters
        ----------
        filename : str
            CSV file name from which to read.
        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : BreedingValueMatrix
            A BreedingValueMatrix read from a CSV file.
        """
        raise NotImplementedError("class method is abstract")
    
    @classmethod
    @abstractmethod
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str]
        ) -> 'BreedingValueMatrix':
        """
        Read a BreedingValueMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which object data is stored.
            If None, object is read from base HDF5 group.

        Returns
        -------
        out : BreedingValueMatrix
            A BreedingValueMatrix read from an HDF5 file.
        """
        raise NotImplementedError("class method is abstract")



################################## Utilities ###################################
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
        raise TypeError("variable '{0}' must be of type '{1}' but received type '{2}'".format(vname,BreedingValueMatrix.__name__,type(v).__name__))
