"""
Module defining basal matrix interfaces and associated error checking routines
for haplotype matrices.
"""

import numbers
from typing import Any, Optional
import numpy
from numpy.typing import DTypeLike
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix

class HaplotypeMatrix(TaxaVariantMatrix,HDF5InputOutput):
    """
    An abstract class for haplotype matrix objects.

    The purpose of this abstract class is to proved base functionality for:
        1) haplotype metadata manipulation routines.
        2) haplotype math functions.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            **kwargs: dict
        ) -> None:
        """
        HaplotypeMatrix constructor

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(HaplotypeMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## General matrix properties ###############
    def ploidy():
        doc = "Ploidy number represented by matrix property."
        def fget(self):
            """Get matrix ploidy number"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set matrix ploidy number"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete matrix ploidy number"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ploidy = property(**ploidy())

    # REMARK: this property is defined in PhasedMatrix as well. the purpose of
    # adding this here as well is to facilitate determination of the number of
    # phases represented by the HaplotypeMatrix. If 0, unphased; if >0, phased.
    def nphase():
        doc = "Number of chromosome phases represented by the matrix."
        def fget(self):
            """Get number of phases"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of phases"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of phases"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    nphase = property(**nphase())

    def mat_format():
        doc = "Matrix representation format property."
        def fget(self):
            """Get matrix representation format"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set matrix representation format"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete matrix representation format"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    mat_format = property(**mat_format())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############## Matrix summary statistics ###############
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

    def meh(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numbers.Number:
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



################################################################################
################################## Utilities ###################################
################################################################################
def is_HaplotypeMatrix(v: Any) -> bool:
    """Return whether an object is a HaplotypeMatrix or not"""
    return isinstance(v, HaplotypeMatrix)

def check_is_HaplotypeMatrix(v: Any, varname: str) -> None:
    """Raise TypeError if object is not a HaplotypeMatrix"""
    if not isinstance(v, HaplotypeMatrix):
        raise TypeError("'%s' must be a HaplotypeMatrix." % varname)
