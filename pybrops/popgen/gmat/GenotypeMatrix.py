"""
Module defining basal matrix interfaces and associated error checking routines
for genotype matrices.
"""

from abc import ABCMeta
from abc import abstractmethod
from numbers import Real
import numpy
from numpy.typing import DTypeLike
from typing import Optional
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybrops.popgen.gmap.GeneticMappableMatrix import GeneticMappableMatrix

class GenotypeMatrix(
        TaxaVariantMatrix,
        GeneticMappableMatrix,
        HDF5InputOutput,
        metaclass = ABCMeta,
    ):
    """
    An abstract class for genoypte matrix objects.

    The purpose of this abstract class is to define base functionality for:
        1) Genotype matrix ploidy and phase metadata.
        2) Genotype matrix format conversion.
        3) Genotype matrix allele counting routines.
        4) Genotype matrix genotype counting routines.
    """

    ########################## Special Object Methods ##########################

    ############################ Object Properties #############################

    ############## Matrix Metadata Properites ##############
    @property
    @abstractmethod
    def ploidy(self) -> int:
        """The ploidy level represented by the genotype matrix."""
        raise NotImplementedError("property is abstract")
    @ploidy.setter
    @abstractmethod
    def ploidy(self, value: int) -> None:
        """Set ploidy level"""
        raise NotImplementedError("property is abstract")
    
    @property
    @abstractmethod
    def nphase(self) -> int:
        """The number of phases represented by the genotype matrix."""
        raise NotImplementedError("property is abstract")
    @nphase.setter
    @abstractmethod
    def nphase(self, value: int) -> None:
        """Set number of phases represented by the genotype matrix"""
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

    ################## Matrix conversion ###################
    @abstractmethod
    def mat_asformat(
            self, 
            format: str
        ) -> numpy.ndarray:
        """
        Get the genotype matrix in a specific format type.

        Parameters
        ----------
        format : str
            Desired output format. Options are ``"{0,1,2}"``, ``"{-1,0,1}"``, ``"{-1,m,1}"``.

        Returns
        -------
        out : numpy.ndarray
            Matrix in the desired output format.
        """
        raise NotImplementedError("method is abstract")

    ############## Matrix summary statistics ###############
    @abstractmethod
    def tacount(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Allele count of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the accumulator and returned array.
            If ``None``, use the native accumulator type (int or float).

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(n,p)`` containing allele counts of the
            allele coded as ``1`` for all ``n`` individuals, for all ``p`` loci.

            Where:

            - ``n`` is the number of taxa (individuals).
            - ``p`` is the number of variants (loci).
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def tafreq(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Allele frequency of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(n,p)`` containing allele frequencies of
            the allele coded as ``1`` for all ``n`` individuals, for all ``p``
            loci.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def acount(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Allele count of the non-zero allele across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele counts of the
            allele coded as ``1`` for all ``p`` loci.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def afreq(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Allele frequency of the non-zero allele across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies of
            the allele coded as ``1`` for all ``p`` loci.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def afixed(
            self,
            dtype: Optional[DTypeLike],
        ) -> numpy.ndarray:
        """
        Determine allele fixation for loci across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.
        
        Returns
        -------
        out : numpy.ndarray
            A ``numpy.ndarray`` of shape ``(p,)`` containing indicator variables 
            for whether a locus is fixed at a particular locus.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def apoly(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Allele polymorphism presence or absense across all loci.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing indicator variables for
            whether the locus is polymorphic.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def maf(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Minor allele frequency across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies for
            the minor allele.
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
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : Real
            A number representing the mean expected heterozygous.
            If ``dtype`` is ``None``, then a native 64-bit floating point is
            returned. Otherwise, of type specified by ``dtype``.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gtcount(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Gather genotype counts for homozygous major, heterozygous, homozygous
        minor for all individuals.

        Parameters
        ----------
        dtype : dtype, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray array of shape ``(g,p)`` containing allele counts
            across all ``p`` loci for each of ``g`` genotype combinations.

            Where:

            - ``out[0]`` is the count of ``0`` genotype across all loci
            - ``out[1]`` is the count of ``1`` genotype across all loci
            - ``out[2]`` is the count of ``2`` genotype across all loci
            - ``...``
            - ``out[g-1]`` is the count of ``g-1`` genotype across all loci
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gtfreq(
            self, 
            dtype: Optional[DTypeLike]
        ) -> numpy.ndarray:
        """
        Gather genotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Parameters
        ----------
        dtype : DTypeLike, None
            The data type of the returned array. 
            If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray array of shape ``(g,p)`` containing haplotype counts
            across all ``p`` loci for all ``g`` genotype combinations.

            Where:

            - ``out[0]`` is the frequency of ``0`` genotype across all loci
            - ``out[1]`` is the frequency of ``1`` genotype across all loci
            - ``out[2]`` is the frequency of ``2`` genotype across all loci
            - ``...``
            - ``out[g-1]`` is the frequency of ``g-1`` genotype across all loci
        """
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GenotypeMatrix(v: object, vname: str) -> None:
    """
    Check if object is of type ``GenotypeMatrix``. Otherwise raise ``TypeError``.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in ``TypeError`` message.
    """
    if not isinstance(v, GenotypeMatrix):
        raise TypeError("'%s' must be a GenotypeMatrix." % vname)

def check_GenotypeMatrix_has_taxa(v: GenotypeMatrix, vname: str) -> None:
    """
    Check if a GenotypeMatrix has taxa labels. Otherwise raise ValueError.

    Parameters
    ----------
    v : GenotypeMatrix
        A GenotypeMatrix object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if v.taxa is None:
        raise ValueError("GenotypeMatrix '{0}' must have a non-None 'taxa' field")
