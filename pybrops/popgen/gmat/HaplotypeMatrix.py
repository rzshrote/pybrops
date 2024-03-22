"""
Module defining basal matrix interfaces and associated error checking routines
for haplotype matrices.
"""

__all__ = [
    "HaplotypeMatrix",
    "check_is_HaplotypeMatrix",
]

from abc import ABCMeta
from abc import abstractmethod
from numbers import Real
from typing import Optional
import numpy
from numpy.typing import DTypeLike
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix

class HaplotypeMatrix(TaxaVariantMatrix,HDF5InputOutput,metaclass=ABCMeta):
    """
    An abstract class for haplotype matrix objects.

    The purpose of this abstract class is to proved base functionality for:
        1) haplotype metadata manipulation routines.
        2) haplotype math functions.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## General matrix properties ###############
    ############## Matrix Metadata Properites ##############
    @property
    @abstractmethod
    def ploidy(self) -> int:
        """The ploidy level represented by the haplotype matrix."""
        raise NotImplementedError("property is abstract")
    @ploidy.setter
    @abstractmethod
    def ploidy(self, value: int) -> None:
        """Set ploidy level"""
        raise NotImplementedError("property is abstract")
    
    # REMARK: this property is defined in PhasedMatrix as well. the purpose of
    # adding this here as well is to facilitate determination of the number of
    # phases represented by the HaplotypeMatrix. If 0, unphased; if >0, phased.
    @property
    @abstractmethod
    def nphase(self) -> int:
        """The number of phases represented by the haplotype matrix."""
        raise NotImplementedError("property is abstract")
    @nphase.setter
    @abstractmethod
    def nphase(self, value: int) -> None:
        """Set number of phases represented by the haplotype matrix"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def mat_format(self) -> str:
        """Matrix representation format."""
        raise NotImplementedError("property is abstract")
    @mat_format.setter
    @abstractmethod
    def mat_format(self, value: str) -> None:
        """Set matrix representation format"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ############## Matrix summary statistics ###############
    @abstractmethod
    def thcount(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Haplotype count of the non-zero haplotype within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, h) containing haplotype counts of the
            haplotype coded as 1 for all 'n' individuals, for all 'h'
            haplotypes.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def thfreq(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Haplotype frequency of the non-zero haplotype within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array and of the accumulator in which the
            elements are summed.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, h) containing haplotype frequencies of
            the haplotype coded as 1 for all 'n' individuals, for all 'h'
            haplotypes.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def hcount(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Haplotype count of the non-zero haplotype across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype counts of the
            haplotype coded as 1 for all 'h' haplotypes.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def hfreq(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Haplotype frequency of the non-zero haplotype across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype frequencies of
            the haplotype coded as 1 for all 'h' haplotypes.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def mhf(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Minor haplotype frequency across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype frequencies for
            the minor haplotype.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def meh(
            self, 
            dtype: Optional[DTypeLike]
        ) -> Real:
        """
        Mean expected heterozygosity across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        meh : numpy.float64
            A 64-bit floating point representing the mean expected
            heterozygosity.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gtcount(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Gather haplotype counts across for homozygous major, heterozygous,
        homozygous minor all individuals.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            An int64 array of shape ``(3,h)`` containing haplotype counts across
            all ``h`` haplotypes.

            Where:

            - ``out[0]`` is the count of ``0`` genotype across all loci
            - ``out[1]`` is the count of ``1`` genotype across all loci
            - ``out[2]`` is the count of ``2`` genotype across all loci
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gtfreq(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Gather haplotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            An float64 array of shape ``(3,h)`` containing haplotype counts
            across all ``h`` haplotypes.

            Where:

            - ``out[0]`` is the frequency of ``0`` genotype across all loci
            - ``out[1]`` is the frequency of ``1`` genotype across all loci
            - ``out[2]`` is the frequency of ``2`` genotype across all loci
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_HaplotypeMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``HaplotypeMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, HaplotypeMatrix):
        raise TypeError("'%s' must be a HaplotypeMatrix." % vname)
