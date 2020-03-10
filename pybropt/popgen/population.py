import cyvcf2
import numpy
import os
import sys
from util.error_subroutines import *
from util.subroutines import *


class Population:
    """docstring for Population."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################

    def __init__(self, geno, genomic_model = None, genetic_map = None,
            taxa = None):
        """
        Population object constructor.

        Parameters
        ----------
        geno : numpy.ndarray
            A int8 binary genotype matrix of shape (m, n, p).
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
        genomic_model : ParametricGenomicModel
            A ParametricGenomicModel object with marker estimates for each
            marker.
        genetic_map : GeneticMap
            A GeneticMap object with positions for each marker.
        taxa : numpy.ndarray, None
            A taxa name matrix of shape (n,).
            Where:
                'n' is the number of individuals.
        """
        # check input data types
        check_is_matrix(geno, "geno")
        cond_check_is_GenomicModel(genomic_model, "genomic_model")
        cond_check_is_GeneticMap(genetic_map, "genetic_map")
        cond_check_is_matrix(taxa, "taxa")

        # check data types
        check_matrix_dtype(geno, "geno", 'int8')
        cond_check_matrix_dtype_is_string_(taxa, "taxa")

        # check matrix number of dimensions
        check_matrix_ndim(geno, "geno", 3)
        cond_check_matrix_ndim(taxa, "taxa", 1)

        # check matrix dimension size alignment
        if taxa is not None:
            check_matrix_axis_len(geno, "geno", 1, len(taxa))

        # set private variables
        self._geno = geno
        self._genomic_model = genomic_model
        self._genetic_map = genetic_map
        self._taxa = taxa

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    def geno():
        doc = "The geno property."
        def fget(self):
            return self._geno
        def fset(self, value):
            self._geno = value
        def fdel(self):
            del self._geno
        return locals()
    geno = property(**geno())

    def taxa():
        doc = "The taxa property."
        def fget(self):
            return self._taxa
        def fset(self, value):
            self._taxa = value
        def fdel(self):
            del self._taxa
        return locals()
    taxa = property(**taxa())

    def genetic_map():
        doc = "The genetic_map property."
        def fget(self):
            return self._genetic_map
        def fset(self, value):
            self._genetic_map = value
        def fdel(self):
            del self._genetic_map
        return locals()
    genetic_map = property(**genetic_map())

    def genomic_model():
        doc = "The genomic_model property."
        def fget(self):
            return self._genomic_model
        def fset(self, value):
            self._genomic_model = value
        def fdel(self):
            del self._genomic_model
        return locals()
    genomic_model = property(**genomic_model())

    def nphase():
        doc = "The nphase property."
        def fget(self):
            return self._geno.shape[0]
        def fset(self, value):
            error_readonly("nphase")
        def fdel(self):
            error_readonly("nphase")
        return locals()
    nphase = property(**nphase())

    def ntaxa():
        doc = "The ntaxa property."
        def fget(self):
            return self._geno.shape[1]
        def fset(self, value):
            error_readonly("ntaxa")
        def fdel(self):
            error_readonly("ntaxa")
        return locals()
    ntaxa = property(**ntaxa())

    def nloci():
        doc = "The nloci property."
        def fget(self):
            return self._geno.shape[2]
        def fset(self, value):
            error_readonly("nloci")
        def fdel(self):
            error_readonly("nloci")
        return locals()
    nloci = property(**nloci())

    ############################################################################
    ############################## Class Methods ###############################
    ############################################################################

    @classmethod
    def D(self):
        raise NotImplementedError("'D' method is not implemented.")

    @classmethod
    def D_gen(self):
        raise NotImplementedError("'D_gen' method is not implemented.")

    @classmethod
    def D_prime(self):
        raise NotImplementedError("'D_prime' method is not implemented.")

    @classmethod
    def D_prime_gen(self):
        raise NotImplementedError("'D_prime_gen' method is not implemented.")

    @classmethod
    def r(self):
        raise NotImplementedError("'r' method is not implemented.")

    @classmethod
    def r_gen(self):
        raise NotImplementedError("'r_gen' method is not implemented.")

    @classmethod
    def r_sq(self):
        raise NotImplementedError("'r_sq' method is not implemented.")

    @classmethod
    def r_sq_gen(self):
        raise NotImplementedError("'r_sq_gen' method is not implemented.")

    @classmethod
    def afreq(self):
        raise NotImplementedError("'afreq' method is not implemented.")

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def from_vcf(fname, base_GenomicModel = None, base_GeneticMap = None,
        kind = 'linear', fill_value = 'extrapolate'):
        # make VCF iterator
        vcf = cyvcf2.VCF(fname)

        # extract taxa names from vcf header
        taxa = numpy.string_(vcf.samples)

        # make empty lists to store extracted values
        geno = []
        chr_grp = []
        chr_start = []
        chr_stop = []
        mkr_name = []

        # iterate through VCF file and accumulate variants
        for variant in vcf:
            # append chromosome string
            chr_grp.append(str(variant.CHROM))

            # append variant position coordinates
            chr_start.append(variant.POS)
            chr_stop.append(variant.POS + len(variant.REF) - 1)

            # append marker name
            mkr_name.append(variant.ID)

            # extract allele states + whether they are phased or not
            phases = numpy.int8(variant.genotypes)

            # check that they are all phased
            check_matrix_all_value(phases[:,2], "is_phased", True)

            # TODO: maybe modify shapes here to avoid transpose and copy below?
            # append genotype states
            geno.append(phases[:,0:2].copy())

        # convert and transpose genotype matrix
        geno = numpy.int8(geno).transpose(2,1,0).copy()

        # convert to object array
        chr_grp = numpy.string_(chr_grp)

        # convert chr_start, chr_stop to numpy.ndarray
        chr_start = numpy.int64(chr_start)
        chr_stop = numpy.int64(chr_stop)

        # make output object
        outpop = Population.from_array(
            geno = geno,
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            taxa = taxa,
            chr_grp = chr_grp,
            chr_pos = chr_pos
        )

        return outpop

    @staticmethod
    def from_array(geno, chr_grp, chr_start, chr_stop, mkr_name = None,
            base_GenomicModel = None, base_GeneticMap = None, taxa = None):
        """
        Align genotypes given several arrays.

        Parameters
        ----------
        geno : numpy.ndarray
        chr_grp : numpy.ndarray
        chr_stop : numpy.ndarray
        mkr_name : numpy.ndarray
        base_GenomicModel : GenomicModel
        base_GeneticMap : GeneticMap

        Returns
        -------
        population : Population
            A Population object.
        """
        # check data types
        check_is_matrix(geno, "geno")
        check_is_matrix(chr_grp, "chr_grp")
        check_is_matrix(chr_start, "chr_start")
        check_is_matrix(chr_stop, "chr_stop")
        cond_check_is_matrix(mkr_name, "mkr_name")
        cond_check_is_GenomicModel(base_GenomicModel, "base_GenomicModel")
        cond_check_is_GeneticMap(base_GeneticMap, "base_GeneticMap")

        # check dimensions
        check_matrix_ndim(geno, "geno", 3)
        check_matrix_ndim(chr_grp, "chr_grp", 1)
        check_matrix_ndim(chr_start, "chr_start", 1)
        check_matrix_ndim(chr_stop, "chr_stop", 1)
        cond_check_matrix_ndim(mkr_name, "mkr_name", 1)

        # check matrix compatiblity lengths
        nloci = geno.shape[2]
        check_matrix_size(chr_grp, "chr_grp", nloci)
        check_matrix_size(chr_start, "chr_start", nloci)
        check_matrix_size(chr_stop, "chr_stop", nloci)
        cond_check_matrix_size(mkr_name, "mkr_name", nloci)

        # check if chr_grp is sorted
        if matrix_is_sorted(chr_grp):
            pass
        else:
            if isinstance(type(base_GenomicModel), NonparametricGenomicModel):
                raise ValueError(
                    "Cannot sort marker data accordingly because of the non-"\
                    "parametric nature of 'base_GenomicModel' "\
                    "(NonparametricGenomicModel does not have a 'coeff' field "\
                    "that can be sorted)."
                )

    @staticmethod
    def preprocess_geno_genetic_map
