# append paths
import sys
import os
popgen_dir = os.path.dirname(os.path.realpath(__file__))    # get pybropt/test/popgen
test_dir = os.path.dirname(popgen_dir)                      # get pybropt/test
pybropt_dir = os.path.dirname(test_dir)                     # get pybropt
sys.path.append(pybropt_dir)                                # append pybropt

# import 3rd party modules we'll need
import cyvcf2
import numpy

# import our libraries
import popgen
import util

class Population:
    """
    Object that represents a breeding population.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, geno, marker_set = None, genomic_model = None, taxa = None):
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
        marker_set : MarkerSet
            A MarkerSet object with positions for each marker.
            This can either be an object of the base class, or any of its
            derived classes (e.g. GeneticMap).
        genomic_model : GenomicModel
            A GenomicModel object for the marker_set.
        taxa : numpy.ndarray, None
            A taxa name matrix of shape (n,).
            Where:
                'n' is the number of individuals.
        """
        # check input data types
        util.check_is_matrix(geno, "geno")
        util.cond_check_is_MarkerSet(marker_set, "marker_set")
        util.cond_check_is_GenomicModel(genomic_model, "genomic_model")
        util.cond_check_is_matrix(taxa, "taxa")

        # check data types
        util.check_matrix_dtype(geno, "geno", 'int8')
        util.cond_check_matrix_dtype_is_string_(taxa, "taxa")

        # check matrix number of dimensions
        util.check_matrix_ndim(geno, "geno", 3)
        util.cond_check_matrix_ndim(taxa, "taxa", 1)

        # check matrix dimension size alignment
        if taxa is not None:
            util.check_matrix_axis_len(geno, "geno", 1, len(taxa))

        # set private variables
        self._geno = geno
        self._genomic_model = genomic_model
        self._marker_set = marker_set
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

    def marker_set():
        doc = "The marker_set property."
        def fget(self):
            return self._marker_set
        def fset(self, value):
            self._marker_set = value
        def fdel(self):
            del self._marker_set
        return locals()
    marker_set = property(**marker_set())

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
    ############################## Object Methods ##############################
    ############################################################################
    def calc_D(self):
        raise NotImplementedError("'D' method is not implemented.")

    def calc_D_prime(self):
        raise NotImplementedError("'D_prime' method is not implemented.")

    def calc_r(self):
        raise NotImplementedError("'r' method is not implemented.")

    def calc_r_sq(self):
        raise NotImplementedError("'r_sq' method is not implemented.")

    def D_gen(self):
        raise NotImplementedError("'D_gen' method is not implemented.")

    def D_prime_gen(self):
        raise NotImplementedError("'D_prime_gen' method is not implemented.")

    def r_gen(self):
        raise NotImplementedError("'r_gen' method is not implemented.")

    def r_sq_gen(self):
        raise NotImplementedError("'r_sq_gen' method is not implemented.")

    def afreq(self, sel = None):
        """
        Calculate allele frequencies for a given subset of the population.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.

        Returns
        -------
        allele_freq : numpy.ndarray
            An allele frequency matrix of shape (p,)
            Where:
                'p' is the number of marker loci.
        """
        # if no selections have been made, select all
        if sel is None:
            sel = slice(None)

        # get selected genotype view
        sgeno = self._geno[:,sel,:]

        # calculate number of chromosome phases in selection
        phases = numpy.float64(sgeno.shape[0] * sgeno.shape[1])

        # calculate allele frequency
        allele_freq = sgeno.sum(0,1) / phases

        return allele_freq

    def gebv(self, sel = None, objcoeff = None):
        """
        Calculate genomic estimated breeding values (GEBVs) using the internal
        genomic_model.

        Parameters
        ----------
        sel : numpy.ndarray
            A selection indices matrix of shape (k,)
            Where:
                'k' is the number of individuals to select.
            Each index indicates which individuals to select.
            Each index in 'sel' represents a single individual's row.
            If 'sel' is None, use all individuals.
        objcoeff : numpy.ndarray, None
            An objective coefficients matrix of shape (t,).
            Where:
                't' is the number of objectives.
            These are used to weigh objectives in the weight sum method.
            If None, do not multiply scores by a weight sum vector.

        Returns
        -------
        gebv : numpy.ndarray
            A gebv matrix of shape (k,) or (k, t).
            Where:
                'k' is the number of individuals to select.
                't' is the number of traits.
        """
        # make sure we have a genomic_model
        if self._genomic_model is None:
            raise RuntimeError("Population has not been provided a GenomicModel.")

        # if no individuals have been selected, select all
        if sel is None:
            sel = slice(None)

        # get view of genotype matrix that is the selections
        sgeno = geno[:,sel,:]

        # calculate GEBVs
        gebv = self._genomic_model.predict(sgeno)

        # take the dot product if necessary
        if objcoeff is not None:
            gebv = gebv.dot(objcoeff)

        return gebv

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def from_array(geno, chr_grp, chr_start, chr_stop, mkr_name = None,
            genomic_model = None, base_GeneticMap = None, taxa = None):
        """
        Construct a Population object from arrays.

        Assumes that GenomicModel is ordered.

        Parameters
        ----------
        geno : numpy.ndarray
        chr_grp : numpy.ndarray
        chr_stop : numpy.ndarray
        mkr_name : numpy.ndarray
        genomic_model : GenomicModel
        base_GeneticMap : GeneticMap
        taxa : numpy.ndarray

        Returns
        -------
        population : Population
            A Population object.
        """
        # check data types
        util.check_is_matrix(geno, "geno")
        util.check_is_matrix(chr_grp, "chr_grp")
        util.check_is_matrix(chr_start, "chr_start")
        util.check_is_matrix(chr_stop, "chr_stop")
        util.cond_check_is_matrix(mkr_name, "mkr_name")
        util.cond_check_is_GeneticMap(base_GeneticMap, "base_GeneticMap")
        util.cond_check_is_GenomicModel(genomic_model, "genomic_model")
        util.cond_check_matrix_dtype_is_string_(taxa, "taxa")

        # check dimensions
        util.check_matrix_ndim(geno, "geno", 3)
        util.check_matrix_ndim(chr_grp, "chr_grp", 1)
        util.check_matrix_ndim(chr_start, "chr_start", 1)
        util.check_matrix_ndim(chr_stop, "chr_stop", 1)
        util.cond_check_matrix_ndim(mkr_name, "mkr_name", 1)
        util.cond_check_matrix_ndim(taxa, "taxa", 1)

        # check matrix compatiblity lengths
        nloci = geno.shape[2]
        util.check_matrix_size(chr_grp, "chr_grp", nloci)
        util.check_matrix_size(chr_start, "chr_start", nloci)
        util.check_matrix_size(chr_stop, "chr_stop", nloci)
        util.cond_check_matrix_size(mkr_name, "mkr_name", nloci)
        util.cond_check_matrix_size(taxa, "taxa", geno.shape[1])


        # TODO: interpolate
        # create genomic models
        marker_set = None
        if base_GeneticMap is not None:
            marker_set = base_GeneticMap.interpolate(
                chr_grp = chr_grp,
                chr_start = chr_start,
                chr_stop = chr_stop,
                mkr_name = mkr_name,
                map_fncode = None,
                kind = 'linear',
                fill_value = 'extrapolate'
            )

        # make population
        population = Population(
            geno = geno,
            marker_set = marker_set,
            genomic_model = genomic_model,
            taxa = taxa
        )

        return population

        # # check if chr_grp is sorted
        # if matrix_is_sorted(chr_grp):
        #     pass
        # else:
        #     if isinstance(type(genomic_model), NonparametricGenomicModel):
        #         raise ValueError(
        #             "Cannot sort marker data accordingly because of the non-"\
        #             "parametric nature of 'genomic_model' "\
        #             "(NonparametricGenomicModel does not have a 'coeff' field "\
        #             "that can be sorted)."
        #         )

    @staticmethod
    def from_vcf(fname, base_GeneticMap = None, genomic_model = None,
        auto_sort = True, auto_mkr_rename = False,
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
            util.check_matrix_all_value(phases[:,2], "is_phased", True)

            # TODO: maybe modify shapes here to avoid transpose and copy below?
            # append genotype states
            geno.append(phases[:,0:2].copy())

        # convert and transpose genotype matrix
        geno = numpy.int8(geno).transpose(2,1,0) # may want to copy()?

        # convert to numpy.ndarray
        chr_grp = numpy.string_(chr_grp)    # convert to string array
        chr_start = numpy.int64(chr_start)  # convert to int64 array
        chr_stop = numpy.int64(chr_stop)    # convert to int64 array
        mkr_name = numpy.string_(mkr_name)  # convert to string array

        # declare output variable
        population = None

        # if base_GeneticMap is None, we construct a MarkerSet object
        if base_GeneticMap is None:
            # make marker set
            marker_set = popgen.MarkerSet(
                chr_grp = chr_grp,
                chr_start = chr_start,
                chr_stop = chr_stop,
                mkr_name = mkr_name,
                auto_sort = auto_sort,
                auto_mkr_rename = auto_mkr_rename
            )

            # make population
            population = Population(
                geno = geno,
                marker_set = marker_set,
                genomic_model = genomic_model,
                taxa = taxa,
            )
        else:
            # make output object
            population = Population.from_array(
                geno = geno,
                chr_grp = chr_grp,
                chr_start = chr_start,
                chr_stop = chr_stop,
                mkr_name = mkr_name,
                genomic_model = genomic_model,
                base_GeneticMap = base_GeneticMap,
                taxa = taxa
            )

        return population
