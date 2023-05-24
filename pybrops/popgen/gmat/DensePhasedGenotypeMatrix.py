"""
Module providing implementations of dense phased genotype matrices and
associated error checking routines.
"""

import copy
import numbers
from typing import Any, Optional
import cyvcf2
import numpy
from numpy.typing import DTypeLike

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_int8
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.mat.DensePhasedTaxaVariantMatrix import DensePhasedTaxaVariantMatrix
from pybrops.popgen.gmat.DenseGenotypeMatrix import DenseGenotypeMatrix
from pybrops.popgen.gmat.PhasedGenotypeMatrix import PhasedGenotypeMatrix

class DensePhasedGenotypeMatrix(DenseGenotypeMatrix,DensePhasedTaxaVariantMatrix,PhasedGenotypeMatrix):
    """
    A concrete class for phased genoypte matrix objects.

    The purpose of this concrete class is to implement functionality for:
        1) Genotype matrix ploidy and phase metadata.
        2) Genotype matrix format conversion.
        3) Genotype matrix allele counting routines.
        4) Genotype matrix genotype counting routines.
        5) Loading phased genotype matrices from VCF and HDF5.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self, 
            mat: numpy.ndarray, 
            taxa: Optional[numpy.ndarray] = None, 
            taxa_grp: Optional[numpy.ndarray] = None, 
            vrnt_chrgrp: Optional[numpy.ndarray] = None, 
            vrnt_phypos: Optional[numpy.ndarray] = None, 
            vrnt_name: Optional[numpy.ndarray] = None, 
            vrnt_genpos: Optional[numpy.ndarray] = None, 
            vrnt_xoprob: Optional[numpy.ndarray] = None, 
            vrnt_hapgrp: Optional[numpy.ndarray] = None, 
            vrnt_hapalt: Optional[numpy.ndarray] = None, 
            vrnt_hapref: Optional[numpy.ndarray] = None, 
            vrnt_mask: Optional[numpy.ndarray] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for the class DensePhasedGenotypeMatrix.

        Parameters
        ----------
        mat : numpy.ndarray
            An int8 haplotype matrix. Must be {0,1,2} format.
        taxa : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa names.
            If ``None``, do not store any taxa name information.
        taxa_grp : numpy.ndarray, None
            A numpy.ndarray of shape ``(n,)`` containing taxa groupings.
            If ``None``, do not store any taxa group information.
        vrnt_chrgrp : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            group labels. If ``None``, do not store any variant chromosome group
            label information.
        vrnt_phypos : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            physical positions. If ``None``, do not store any variant chromosome
            physical position information.
        vrnt_name : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant names.
            If ``None``, do not store any variant names.
        vrnt_genpos : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant chromosome
            genetic positions. If ``None``, do not store any variant chromosome
            genetic position information.
        vrnt_xoprob : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant crossover
            probabilities. If ``None``, do not store any variant crossover
            probabilities.
        vrnt_hapgrp : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant haplotype
            group labels. If ``None``, do not store any variant haplotype group
            label information.
        vrnt_hapalt : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant alternative
            alleles. If ``None``, do not store any variant alternative allele
            information.
        vrnt_hapref : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing variant reference
            alleles. If ``None``, do not store any variant reference allele
            information.
        vrnt_mask : numpy.ndarray, None
            A numpy.ndarray of shape ``(p,)`` containing a variant mask.
            If ``None``, do not store any variant mask information.
        """
        # TODO: Add a mat_format option to store as {0,1,2}, {-1,0,1}, etc.
        super(DensePhasedGenotypeMatrix, self).__init__(
            mat = mat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            vrnt_genpos = vrnt_genpos,
            vrnt_xoprob = vrnt_xoprob,
            vrnt_hapgrp = vrnt_hapgrp,
            vrnt_hapalt = vrnt_hapalt,
            vrnt_hapref = vrnt_hapref,
            vrnt_mask = vrnt_mask,
            **kwargs
        )

    #################### Matrix copying ####################
    def __copy__(
            self
        ) -> 'DensePhasedGenotypeMatrix':
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : DensePhasedGenotypeMatrix
            A copyt of the DensePhasedGenotypeMatrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.copy(self.mat),
            taxa = copy.copy(self.taxa),
            taxa_grp = copy.copy(self.taxa_grp),
            vrnt_chrgrp = copy.copy(self.vrnt_chrgrp),
            vrnt_phypos = copy.copy(self.vrnt_phypos),
            vrnt_name = copy.copy(self.vrnt_name),
            vrnt_genpos = copy.copy(self.vrnt_genpos),
            vrnt_xoprob = copy.copy(self.vrnt_xoprob),
            vrnt_hapgrp = copy.copy(self.vrnt_hapgrp),
            vrnt_hapalt = copy.copy(self.vrnt_hapalt),
            vrnt_hapref = copy.copy(self.vrnt_hapref),
            vrnt_mask = copy.copy(self.vrnt_mask)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.copy(self.taxa_grp_name)
        out.taxa_grp_stix = copy.copy(self.taxa_grp_stix)
        out.taxa_grp_spix = copy.copy(self.taxa_grp_spix)
        out.taxa_grp_len = copy.copy(self.taxa_grp_len)

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.copy(self.vrnt_chrgrp_name)
        out.vrnt_chrgrp_stix = copy.copy(self.vrnt_chrgrp_stix)
        out.vrnt_chrgrp_spix = copy.copy(self.vrnt_chrgrp_spix)
        out.vrnt_chrgrp_len = copy.copy(self.vrnt_chrgrp_len)

        return out

    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'DensePhasedGenotypeMatrix':
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict
            Deep copy metadata.

        Returns
        -------
        out : DensePhasedGenotypeMatrix
            A deep copy of the DensePhasedGenotypeMatrix.
        """
        # create new object
        out = self.__class__(
            mat = copy.deepcopy(self.mat, memo),
            taxa = copy.deepcopy(self.taxa, memo),
            taxa_grp = copy.deepcopy(self.taxa_grp, memo),
            vrnt_chrgrp = copy.deepcopy(self.vrnt_chrgrp, memo),
            vrnt_phypos = copy.deepcopy(self.vrnt_phypos, memo),
            vrnt_name = copy.deepcopy(self.vrnt_name, memo),
            vrnt_genpos = copy.deepcopy(self.vrnt_genpos, memo),
            vrnt_xoprob = copy.deepcopy(self.vrnt_xoprob, memo),
            vrnt_hapgrp = copy.deepcopy(self.vrnt_hapgrp, memo),
            vrnt_hapalt = copy.deepcopy(self.vrnt_hapalt, memo),
            vrnt_hapref = copy.deepcopy(self.vrnt_hapref, memo),
            vrnt_mask = copy.deepcopy(self.vrnt_mask, memo)
        )

        # copy taxa metadata
        out.taxa_grp_name = copy.deepcopy(self.taxa_grp_name, memo)
        out.taxa_grp_stix = copy.deepcopy(self.taxa_grp_stix, memo)
        out.taxa_grp_spix = copy.deepcopy(self.taxa_grp_spix, memo)
        out.taxa_grp_len = copy.deepcopy(self.taxa_grp_len, memo)

        # copy variant metadata
        out.vrnt_chrgrp_name = copy.deepcopy(self.vrnt_chrgrp_name, memo)
        out.vrnt_chrgrp_stix = copy.deepcopy(self.vrnt_chrgrp_stix, memo)
        out.vrnt_chrgrp_spix = copy.deepcopy(self.vrnt_chrgrp_spix, memo)
        out.vrnt_chrgrp_len = copy.deepcopy(self.vrnt_chrgrp_len, memo)

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genotype Data Properites ###############
    @DenseGenotypeMatrix.mat.setter
    def mat(self, value: numpy.ndarray) -> None:
        check_is_ndarray(value, "mat")
        check_ndarray_dtype_is_int8(value, "mat")
        check_ndarray_ndim(value, "mat", 3)
        self._mat = value

    ############## General matrix properties ###############
    @DenseGenotypeMatrix.ploidy.getter
    def ploidy(self) -> int:
        """Get matrix ploidy number"""
        return self._mat.shape[self.phase_axis]

    @DenseGenotypeMatrix.mat_format.getter
    def mat_format(self):
        """Get matrix representation format"""
        return "{0,1,2}"

    ############## Phase Metadata Properites ###############
    # this property must be overwritten to what is in DensePhasedMatrix since
    # DenseGenotypeMatrix overwrites it.
    @DenseGenotypeMatrix.nphase.getter
    def nphase(self):
        """Get number of phases"""
        return self._mat.shape[self.phase_axis]

    @DensePhasedTaxaVariantMatrix.phase_axis.getter
    def phase_axis(self):
        """Get phase axis number"""
        return 0

    ############### Taxa Metadata Properites ###############
    @DenseGenotypeMatrix.taxa_axis.getter
    def taxa_axis(self):
        """Get taxa axis number"""
        return 1

    ############# Variant Metadata Properites ##############
    @DenseGenotypeMatrix.vrnt_axis.getter
    def vrnt_axis(self):
        """Get variant axis"""
        return 2

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################## Matrix conversion ###################
    def mat_asformat(
            self, 
            format: str
        ) -> numpy.ndarray:
        """
        Get mat in a specific format type.

        Parameters
        ---------
        format : str
            Desired output format. Options are "{0,1,2}", "{-1,0,1}", "{-1,m,1}".

        Returns
        -------
        out : numpy.ndarray
            Matrix in the desired output format.
        """
        if format == "{0,1,2}":
            return self.mat.sum(0, dtype = self.mat.dtype)
        elif format == "{-1,0,1}":
            out = self.mat.sum(0, dtype = self.mat.dtype)
            out -= 1
            return out
        elif format == "{-1,m,1}":
            # OPTIMIZE: there's probably a matrix multiplication way to do this instead of if-else
            out = self.mat.sum(0) - 1.0    # (n,p) float64
            for i in range(out.shape[1]):
                view = out[:,i]
                mean = view.mean()
                mask = (view == 0)
                out[mask,i] = mean
            return out
        else:
            raise ValueError('Format not recognized. Options are "{0,1,2}", "{-1,0,1}", "{-1,m,1}".')

    ############## Matrix summary statistics ###############
    def tacount(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele count of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the accumulator and returned array.
            If ``None``, use the native accumulator type (int or float).

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(n,p)`` containing allele counts of the
            allele coded as ``1`` for all ``n`` individuals, for all ``p`` loci.
        """
        # get accumulator type
        if dtype is None:
            if numpy.issubdtype(self._mat.dtype, numpy.integer):
                dtype = int
            elif numpy.issubdtype(self._mat.dtype, numpy.floating):
                dtype = float
            else:
                raise ValueError("No default accumulator type for GenotypeMatrix dtype {0}".format(self._mat.dtype))
        # calculate the taxa allele sums
        out = self._mat.sum(
            axis = self.phase_axis, 
            dtype = dtype
        )
        return out

    def tafreq(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele frequency of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(n,p)`` containing allele frequencies of
            the allele coded as ``1`` for all ``n`` individuals, for all ``p``
            loci.
        """
        rnphase = 1.0 / self.ploidy                     # take the reciprocal of the ploidy number
        out = rnphase * self._mat.sum(self.phase_axis)  # take sum across the phase axis (0) and divide by ploidy
        if dtype is not None:                           # if dtype is specified
            dtype = numpy.dtype(dtype)                  # ensure conversion to dtype class
            if out.dtype != dtype:                      # if output dtype and desired are different
                out = dtype.type(out)                   # convert to correct dtype
        return out

    def acount(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele count of the non-zero allele across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele counts of the
            allele coded as ``1`` for all ``p`` loci.
        """
        out = self._mat.sum((self.phase_axis,self.taxa_axis))   # take sum across the phase and taxa axis
        if dtype is not None:                                   # if dtype is specified
            dtype = numpy.dtype(dtype)                          # ensure conversion to dtype class
            if out.dtype != dtype:                              # if output dtype and desired are different
                out = dtype.type(out)                           # convert to correct dtype
        return out

    def afreq(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele frequency of the non-zero allele across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies of
            the allele coded as ``1`` for all ``p`` loci.
        """
        rnphase = 1.0 / (self.ploidy * self.ntaxa)                      # take 1 / (ploidy * ntaxa)
        out = rnphase * self._mat.sum((self.phase_axis,self.taxa_axis)) # take sum across the taxa axis (0) and divide by nphase
        if dtype is not None:                                           # if dtype is specified
            dtype = numpy.dtype(dtype)                                  # ensure conversion to dtype class
            if out.dtype != dtype:                                      # if output dtype and desired are different
                out = dtype.type(out)                                   # convert to correct dtype
        return out

    def apoly(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Allele polymorphism presence or absense across all loci.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing indicator variables for
            whether the locus is polymorphic.
        """
        # find any loci that are monomorphic for 0 allele
        min_mask = numpy.all(self.mat == 0, axis = (self.phase_axis,self.taxa_axis))

        # find any loci that are monomorphic for 1 allele
        max_mask = numpy.all(self.mat == 1, axis = (self.phase_axis,self.taxa_axis))

        # logical or results and take logical not
        out = numpy.logical_not(min_mask | max_mask)

        # convert to specific dtype if needed
        if dtype is not None:                           # if dtype is specified
            dtype = numpy.dtype(dtype)                  # ensure conversion to dtype class
            if out.dtype != dtype:                      # if output dtype and desired are different
                out = dtype.type(out)                   # convert to correct dtype

        return out

    def maf(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Minor allele frequency across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies for
            the minor allele.
        """
        out = self.afreq(dtype)     # get allele frequencies
        mask = out > 0.5            # create mask of allele frequencies > 0.5
        out[mask] = 1.0 - out[mask] # take 1 - allele frequency
        return out

    def meh(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numbers.Number:
        """
        Mean expected heterozygosity across all taxa.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.float64, other
            A number representing the mean expected heterozygous.
            If ``dtype`` is ``None``, then a native 64-bit floating point is
            returned. Otherwise, of type specified by ``dtype``.
        """
        p = self.afreq()                    # get haplotype frequency (p)
        out = (p * (1.0 - p)).sum()         # take p*(1-p) across all loci and sum the products
        out *= (self.ploidy / self.nvrnt)   # multiply summation by (ploidy / nvrnt)
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def gtcount(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Gather genotype counts for homozygous major, heterozygous, homozygous
        minor for all individuals.

        Parameters
        ----------
        dtype : dtype, None
            The dtype of the returned array. If ``None``, use the native type.

        Returns
        -------
        out : numpy.ndarray
            An ``int64`` array of shape ``(g,p)`` containing allele counts across
            all ``p`` loci for each of ``g`` genotype combinations.

            Where:

            - ``out[0]`` is the count of ``0`` genotype across all loci
            - ``out[1]`` is the count of ``1`` genotype across all loci
            - ``out[2]`` is the count of ``2`` genotype across all loci
            - ``...``
            - ``out[g-1]`` is the count of ``g-1`` genotype across all loci
        """
        ngt = self.nphase + 1                   # get number of genotype combos
        mat = self._mat.sum(self.phase_axis)    # get sum of all phases

        out = numpy.empty(      # allocate output array
            (ngt, self.nvrnt),  # shape (g, p)
            dtype='int64'       # int64 dtype
        )

        # get correct axis to sum across
        axis = (self.taxa_axis - 1) if (self.phase_axis < self.taxa_axis) else self.taxa_axis

        for i in range(ngt):
            out[i] = (mat == i).sum(axis)   # record counts for genotype 'i'

        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype

        return out

    def gtfreq(
            self, 
            dtype: Optional[DTypeLike] = None
        ) -> numpy.ndarray:
        """
        Gather genotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Returns
        -------
        out : numpy.ndarray
            An ``float64`` array of shape ``(g,p)`` containing haplotype counts
            across all ``p`` loci for all ``g`` genotype combinations.

            Where:

            - ``out[0]`` is the frequency of ``0`` genotype across all loci
            - ``out[1]`` is the frequency of ``1`` genotype across all loci
            - ``out[2]`` is the frequency of ``2`` genotype across all loci
            - ``...``
            - ``out[g-1]`` is the frequency of ``g-1`` genotype across all loci
        """
        recip = 1.0 / self.ntaxa        # get reciprocal of number of taxa
        out = recip * self.gtcount()    # calculate genotype frequencies
        if dtype is not None:           # if dtype is specified
            dtype = numpy.dtype(dtype)  # ensure conversion to dtype class
            if out.dtype != dtype:      # if output dtype and desired are different
                out = dtype.type(out)   # convert to correct dtype
        return out

    ################### Matrix File I/O ####################

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    ################### Matrix File I/O ####################
    @classmethod
    def from_vcf(
            cls, 
            fname: str
        ) -> 'DensePhasedGenotypeMatrix':
        """
        Does not ensure that data is phased, just reads it as phased.
        """
        # make VCF iterator
        vcf = cyvcf2.VCF(fname)

        # extract taxa names from vcf header
        taxa = numpy.object_(vcf.samples)

        # make empty lists to store extracted values
        mat = []
        vrnt_chrgrp = []
        vrnt_phypos = []
        vrnt_name = []

        # iterate through VCF file and accumulate variants
        for variant in vcf:
            # append chromosome integer
            vrnt_chrgrp.append(int(variant.CHROM))

            # append variant position coordinates
            vrnt_phypos.append(variant.POS)

            # append marker name
            vrnt_name.append(str(variant.ID))

            # extract allele states + whether they are phased or not
            phases = numpy.int8(variant.genotypes)

            # check that they are all phased
            #pybrops.util.check_matrix_all_value(phases[:,2], "is_phased", True)

            # TODO: maybe modify shapes here to avoid transpose and copy below?
            # append genotype states
            mat.append(phases[:,0:2].copy())

        # convert and transpose genotype matrix
        mat = numpy.int8(mat).transpose(2,1,0) # may want to copy()?

        # convert to numpy.ndarray
        vrnt_chrgrp = numpy.int64(vrnt_chrgrp)  # convert to int64 array
        vrnt_phypos = numpy.int64(vrnt_phypos)  # convert to int64 array
        vrnt_name = numpy.object_(vrnt_name)    # convert to object array

        out = cls(
            mat = mat,
            vrnt_chrgrp = vrnt_chrgrp,
            vrnt_phypos = vrnt_phypos,
            vrnt_name = vrnt_name,
            taxa = taxa
        )

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedGenotypeMatrix(v: object) -> bool:
    return isinstance(v, DensePhasedGenotypeMatrix)

def check_is_DensePhasedGenotypeMatrix(v: object, vname: str) -> None:
    if not isinstance(v, DensePhasedGenotypeMatrix):
        raise TypeError("'{0}' must be a DensePhasedGenotypeMatrix.".format(varname))
