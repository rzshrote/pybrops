import copy
import numpy
import cyvcf2

from . import PhasedGenotypeMatrix
from . import DenseGenotypeMatrix
from pybropt.core.mat import DensePhasedTaxaVariantMatrix
from pybropt.core.mat import get_axis

from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_dtype_is_int8
from pybropt.core.error import check_ndarray_is_3d
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_ndarray_dtype_is_object
from pybropt.core.error import cond_check_ndarray_ndim
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import cond_check_ndarray_dtype_is_int64
from pybropt.core.error import cond_check_ndarray_dtype_is_float64
from pybropt.core.error import cond_check_ndarray_dtype_is_bool
from pybropt.core.error import check_is_int
from pybropt.core.error import error_readonly

from pybropt.popgen.gmap import check_is_GeneticMap
from pybropt.popgen.gmap import check_is_GeneticMapFunction

class DensePhasedGenotypeMatrix(DenseGenotypeMatrix,DensePhasedTaxaVariantMatrix,PhasedGenotypeMatrix):
    """docstring for DensePhasedGenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, taxa = None, taxa_grp = None, vrnt_chrgrp = None, vrnt_phypos = None, vrnt_name = None, vrnt_genpos = None, vrnt_xoprob = None, vrnt_hapgrp = None, vrnt_hapalt = None, vrnt_hapref = None, vrnt_mask = None, **kwargs):
        """
        Parameters
        ----------
        mat : numpy.ndarray
            An int8 haplotype matrix. Must be {0,1,2} format.
        ploidy : int
            The ploidy represented by the haplotype matrix. This only represents
            ploidy of the reproductive habit. If the organism represented is an
            allopolyploid (e.g. hexaploid wheat), the ploidy is 2 since it
            reproduces in a diploid manner.
        # TODO: Add a mat_format option to store as {0,1,2}, {-1,0,1}, etc.
        """
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
    def __copy__(self):
        """
        Make a shallow copy of the the matrix.

        Returns
        -------
        out : Matrix
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

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the matrix.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : Matrix
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

    ############## Genotype Data Properites ##############
    def mat():
        doc = "The mat property."
        def fget(self):
            return self._mat
        def fset(self, value):
            check_is_ndarray(value, "mat")
            check_ndarray_dtype_is_int8(value, "mat")
            check_ndarray_is_3d(value, "mat")
            self._mat = value
        def fdel(self):
            del self._mat
        return locals()
    mat = property(**mat())

    ############## General matrix properties ###############
    def ploidy():
        doc = "Ploidy number represented by matrix property."
        def fget(self):
            """Get matrix ploidy number"""
            return self._mat.shape[self.phase_axis]
        def fset(self, value):
            """Set matrix ploidy number"""
            error_readonly("ploidy")
        def fdel(self):
            """Delete matrix ploidy number"""
            error_readonly("ploidy")
        return locals()
    ploidy = property(**ploidy())

    def mat_format():
        doc = "Matrix representation format property."
        def fget(self):
            """Get matrix representation format"""
            return "{0,1,2}"
        def fset(self, value):
            """Set matrix representation format"""
            error_readonly("mat_format")
        def fdel(self):
            """Delete matrix representation format"""
            error_readonly("mat_format")
        return locals()
    mat_format = property(**mat_format())

    ############## Phase Metadata Properites ###############
    # this property must be overwritten to what is in DensePhasedMatrix since
    # DenseGenotypeMatrix overwrites it.
    def nphase():
        doc = "Number of chromosome phases represented by the matrix."
        def fget(self):
            """Get number of phases"""
            return self._mat.shape[self.phase_axis]
        def fset(self, value):
            """Set number of phases"""
            error_readonly("nphase")
        def fdel(self):
            """Delete number of phases"""
            error_readonly("nphase")
        return locals()
    nphase = property(**nphase())

    def phase_axis():
        doc = "Axis along which phases are stored property."
        def fget(self):
            """Get phase axis number"""
            return 0
        def fset(self, value):
            """Set phase axis number"""
            error_readonly("phase_axis")
        def fdel(self):
            """Delete phase axis number"""
            error_readonly("phase_axis")
        return locals()
    phase_axis = property(**phase_axis())

    ############### Taxa Metadata Properites ###############
    def taxa_axis():
        doc = "Axis along which taxa are stored property."
        def fget(self):
            """Get taxa axis number"""
            return 1
        def fset(self, value):
            """Set taxa axis number"""
            error_readonly("taxa_axis")
        def fdel(self):
            """Delete taxa axis number"""
            error_readonly("taxa_axis")
        return locals()
    taxa_axis = property(**taxa_axis())

    ############# Variant Metadata Properites ##############
    def vrnt_axis():
        doc = "Axis along which variants are stored property."
        def fget(self):
            """Get variant axis"""
            return 2
        def fset(self, value):
            """Set variant axis"""
            error_readonly("vrnt_axis")
        def fdel(self):
            """Delete variant axis"""
            error_readonly("vrnt_axis")
        return locals()
    vrnt_axis = property(**vrnt_axis())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################## Matrix conversion ###################
    def mat_asformat(self, format):
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
    def tacount(self, dtype = None):
        """
        Allele count of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array and of the accumulator in which the
            elements are summed.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, p) containing allele counts of the
            allele coded as 1 for all 'n' individuals, for all 'p' loci.
        """
        out = self._mat.sum(self.phase_axis)    # take sum across the phase axis
        if dtype is not None:                   # if dtype is specified
            dtype = numpy.dtype(dtype)          # ensure conversion to dtype class
            if out.dtype != dtype:              # if output dtype and desired are different
                out = dtype.type(out)           # convert to correct dtype
        return out

    def tafreq(self, dtype = None):
        """
        Allele frequency of the non-zero allele within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array and of the accumulator in which the
            elements are summed.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, p) containing allele frequencies of the
            allele coded as 1 for all 'n' individuals, for all 'p' loci.
        """
        rnphase = 1.0 / self.ploidy                     # take the reciprocal of the ploidy number
        out = rnphase * self._mat.sum(self.phase_axis)  # take sum across the phase axis (0) and divide by ploidy
        if dtype is not None:                           # if dtype is specified
            dtype = numpy.dtype(dtype)                  # ensure conversion to dtype class
            if out.dtype != dtype:                      # if output dtype and desired are different
                out = dtype.type(out)                   # convert to correct dtype
        return out

    def acount(self, dtype = None):
        """
        Allele count of the non-zero allele across all taxa.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele counts of the allele
            coded as 1 for all 'p' loci.
        """
        out = self._mat.sum((self.phase_axis,self.taxa_axis))   # take sum across the phase and taxa axis
        if dtype is not None:                                   # if dtype is specified
            dtype = numpy.dtype(dtype)                          # ensure conversion to dtype class
            if out.dtype != dtype:                              # if output dtype and desired are different
                out = dtype.type(out)                           # convert to correct dtype
        return out

    def afreq(self, dtype = None):
        """
        Allele frequency of the non-zero allele across all taxa.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele frequencies of the
            allele coded as 1 for all 'p' loci.
        """
        rnphase = 1.0 / (self.ploidy * self.ntaxa)                      # take 1 / (ploidy * ntaxa)
        out = rnphase * self._mat.sum((self.phase_axis,self.taxa_axis)) # take sum across the taxa axis (0) and divide by nphase
        if dtype is not None:                                           # if dtype is specified
            dtype = numpy.dtype(dtype)                                  # ensure conversion to dtype class
            if out.dtype != dtype:                                      # if output dtype and desired are different
                out = dtype.type(out)                                   # convert to correct dtype
        return out

    def maf(self, dtype = None):
        """
        Minor allele frequency across all taxa.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (p,) containing allele frequencies for the
            minor allele.
        """
        out = self.afreq(dtype)     # get allele frequencies
        mask = out > 0.5            # create mask of allele frequencies > 0.5
        out[mask] = 1.0 - out[mask] # take 1 - allele frequency
        return out

    def mehe(self):
        """
        Mean expected heterozygosity across all taxa.

        Returns
        -------
        mehe : numpy.float64
            A 64-bit floating point representing the mean expected
            heterozygosity.
        """
        p = self.afreq()                    # get haplotype frequency (p)
        out = (p * (1.0 - p)).sum()         # take p*(1-p) across all loci and sum the products
        out *= (self.ploidy / self.nvrnt)   # multiply summation by (nphase / nvrnt)
        return numpy.float64(out)

    def gtcount(self):
        """
        Gather genotype counts for homozygous major, heterozygous, homozygous
        minor for all individuals.

        Returns
        -------
        out : numpy.ndarray
            An int64 array of shape (g, p) containing allele counts across
            all 'p' loci for each of 'g' genotype combinations.
            Rows are as follows:
                out[0] = count of '0' genotype across all loci
                out[1] = count of '1' genotype across all loci
                out[2] = count of '2' genotype across all loci
                ...
                out[g-1] = count of 'g-1' genotype across all loci
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

        return out

    def gtfreq(self):
        """
        Gather genotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Returns
        -------
        out : numpy.ndarray
            An float64 array of shape (g, p) containing haplotype counts across
            all 'p' loci for all 'g' genotype combinations.
            Rows are as follows:
                out[0] = frequency of '0' genotype across all loci
                out[1] = frequency of '1' genotype across all loci
                out[2] = frequency of '2' genotype across all loci
                ...
                out[g-1] = frequency of 'g-1' genotype across all loci
        """
        recip = 1.0 / self.ntaxa        # get reciprocal of number of taxa
        out = recip * self.gtcount()    # calculate genotype frequencies
        return out

    @classmethod
    def from_vcf(cls, fname):
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
            #pybropt.util.check_matrix_all_value(phases[:,2], "is_phased", True)

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
def is_DensePhasedGenotypeMatrix(v):
    return isinstance(v, DensePhasedGenotypeMatrix)

def check_is_DensePhasedGenotypeMatrix(v, varname):
    if not isinstance(v, DensePhasedGenotypeMatrix):
        raise TypeError("'{0}' must be a DensePhasedGenotypeMatrix.".format(varname))

def cond_check_is_DensePhasedGenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DensePhasedGenotypeMatrix(v, varname)
