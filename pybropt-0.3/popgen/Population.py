# import 3rd party modules we'll need
import cyvcf2
import numpy

# import our libraries
import pybropt.popgen.MarkerSet
import pybropt.popgen.GeneticMap
import pybropt.util

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
                Cannot be None.
        marker_set : MarkerSet
            A MarkerSet object with positions for each marker. This can either
            be an object of the base class, or any of its derived classes
            (e.g. GeneticMap).
        genomic_model : GenomicModel
            A GenomicModel object for the marker_set.
        taxa : numpy.ndarray, None
            A taxa name matrix of shape (n,).
            Where:
                'n' is the number of individuals.
            This must be a numpy.string_ matrix.
        """
        # set private variables with error check; this is order dependent!!!
        self.geno = geno                    # geno before all
        self.marker_set = marker_set        # marker_set before genomic_model
        self.genomic_model = genomic_model  # genomic_model last
        self.taxa = taxa

        # set private variables without error check
        self._sorted = False

    def __len__(self):
        """
        Return the number of taxa in the genotype matrix.

        Returns
        -------
        length : int
            Length of taxa axis in the genotype matrix.
        """
        return self._geno.shape[1]

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    def geno():
        doc = "Genotype matrix access. See geno.fget, geno.fset, geno.fdel function docs."
        def fget(self):
            """
            Retrieve a pointer to the genotype matrix.

            Returns
            -------
            geno : numpy.ndarray
                Pointer to the genotype matrix.
            """
            return self._geno
        def fset(self, geno):
            """
            Set the pointer to the genotype matrix. Perform error checks.

            Parameters
            ----------
            geno : numpy.ndarray
                A int8 binary genotype matrix of shape (m, n, p).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'p' is the number of markers.
                Cannot be None.
            """
            # check input data types, then set variable
            pybropt.util.check_is_matrix(geno, "geno")
            pybropt.util.check_matrix_dtype(geno, "geno", 'int8')
            pybropt.util.check_matrix_ndim(geno, "geno", 3)
            self._geno = geno
        def fdel(self):
            """
            Remove the genotype matrix (self._geno) from scope.
            """
            del self._geno
        return locals()
    geno = property(**geno())

    def taxa():
        doc = "Taxa name matrix access. See taxa.fget, taxa.fset, taxa.fdel function docs."
        def fget(self):
            """
            Retrieve a pointer to the taxa name matrix.

            Returns
            -------
            taxa : numpy.ndarray
                Pointer to the taxa name matrix.
            """
            return self._taxa
        def fset(self, taxa):
            """
            Set the pointer to the taxa name matrix. Perform error checks.

            Parameters
            ----------
            taxa : numpy.ndarray, None
                A taxa name matrix of shape (n,).
                Where:
                    'n' is the number of individuals.
                This must be a numpy.string_ matrix.
            """
            # check input data type if not None, then set variable
            if taxa is not None:
                pybropt.util.check_is_matrix(taxa, "taxa")
                pybropt.util.check_matrix_ndim(taxa, "taxa", 1)
                pybropt.util.check_matrix_dtype_is_string_(taxa, "taxa")
                pybropt.util.check_matrix_size(taxa, "taxa", self._geno.shape[1])
            self._taxa = taxa
        def fdel(self):
            """
            Remove the taxa name matrix (self._taxa) from scope.
            """
            del self._taxa
        return locals()
    taxa = property(**taxa())

    def marker_set():
        doc = "Marker set object access. See marker_set.fget, marker_set.fset, marker_set.fdel function docs."
        def fget(self):
            """
            Retrieve a pointer to the MarkerSet object.

            Returns
            -------
            marker_set : MarkerSet
                Pointer to the MarkerSet object
            """
            return self._marker_set
        def fset(self, marker_set):
            """
            Set the pointer to the marker set object. Perform error checks.

            Parameters
            ----------
            marker_set : MarkerSet
                A MarkerSet object with positions for each marker. This can
                either be an object of the base class, or any of its derived
                classes (e.g. GeneticMap).
            """
            # check input data type if not None, then set variable
            if marker_set is not None:
                pybropt.util.check_is_MarkerSet(marker_set, "marker_set")
                pybropt.util.check_matrix_axis_len(self._geno, "self.geno", 2, len(marker_set))
            self._marker_set = marker_set
        def fdel(self):
            """
            Remove the MarkerSet object (self._marker_set) from scope.
            """
            del self._marker_set
        return locals()
    marker_set = property(**marker_set())

    def genomic_model():
        doc = "Genomic model object access. See genomic_model.fget, genomic_model.fset, genomic_model.fdel function docs."
        def fget(self):
            """
            Retrieve a pointer to the GenomicModel object.

            Returns
            -------
            genomic_model : GenomicModel
                Pointer to the GenomicModel object.
            """
            return self._genomic_model
        def fset(self, genomic_model):
            """
            Set the pointer to the marker set object. Perform error checks.

            Parameters
            ----------
            genomic_model : GenomicModel
                A GenomicModel object with positions for each marker. This can
                either be an object of the base class, or any of its derived
                classes (e.g. ParametricGenomicModel).
            """
            # if check data type if not None, and set value
            if genomic_model is not None:
                pybropt.util.check_is_GenomicModel(genomic_model, "genomic_model")
                # if we have a marker set, align the genotypes and markers
                if self._marker_set is not None:
                    indices = genomic_model.mkr_where(self._marker_set.mkr_name)
                    genomic_model.reorder(indices)
            self._genomic_model = genomic_model
        def fdel(self):
            """
            Remove the MarkerSet object (self._marker_set) from scope.
            """
            del self._genomic_model
        return locals()
    genomic_model = property(**genomic_model())

    def nphase():
        doc = "Number of chromosome phases (the ploidy). [readonly]"
        def fget(self):
            """
            Retrieve the number of chromosome phases (the ploidy of the genotype data).

            Returns
            -------
            nphase : int
                The number of chromosome phases according to the genotype matrix shape.
            """
            return self._geno.shape[0]
        def fset(self, value):
            """
            Set the number of chromosome phases.
            Always raises an AttributeError. This property is readonly!
            """
            pybropt.util.error_readonly("nphase")
        def fdel(self):
            """
            Remove the number of chromosome phases from scope.
            Always raises an AttributeError. This property is readonly!
            """
            pybropt.util.error_readonly("nphase")
        return locals()
    nphase = property(**nphase())

    def ntaxa():
        doc = "Number of genotype taxa. [readonly]"
        def fget(self):
            """
            Retrieve the number of genotype taxa.

            Returns
            -------
            ntaxa : int
                Number of genotype taxa according to the genotype matrix shape.
            """
            return self._geno.shape[1]
        def fset(self, value):
            """
            Set the number of genotype taxa.
            Always raises an AttributeError. This property is readonly!
            """
            pybropt.util.error_readonly("ntaxa")
        def fdel(self):
            """
            Remove the number of genotype taxa from scope.
            Always raises an AttributeError. This property is readonly!
            """
            pybropt.util.error_readonly("ntaxa")
        return locals()
    ntaxa = property(**ntaxa())

    def nloci():
        doc = "Number of marker loci. [readonly]"
        def fget(self):
            """
            Retrieve the number of marker loci in the genotype matrix.

            Returns
            -------
            nloci : int
                Number of marker loci according to the genotype matrix shape.
            """
            return self._geno.shape[2]
        def fset(self, value):
            """
            Set the number of marker loci.
            Always raises an AttributeError. This property is readonly!
            """
            pybropt.util.error_readonly("nloci")
        def fdel(self):
            """
            Remove the number of marker loci from scope.
            Always raises an AttributeError. This property is readonly!
            """
            pybropt.util.error_readonly("nloci")
        return locals()
    nloci = property(**nloci())

    def sorted():
        doc = "Boolean variable indicating whether markers have been sorted."
        def fget(self):
            """
            Retrieve whether markers in this object have been sorted.

            Returns
            -------
            sorted : bool
                Whether the Population object has been sorted.
            """
            return self._sorted
        def fset(self, sorted):
            """
            Set whether the Population object has been sorted.
            Performs error checks.

            Parameters
            ----------
            sorted : bool
                Whether the Population object has been sorted.
            """
            # check data types and set variable
            pybropt.util.check_is_bool(sorted, "sorted")
            self._sorted = sorted
        def fdel(self):
            """
            Remove the sorted property (self._sorted) from scope.
            """
            del self._sorted
        return locals()
    sorted = property(**sorted())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def _check_self_has_MarkerSet(self):
        """
        Determine if the Population has a MarkerSet. Raise RuntimeError otherwise.
        """
        if not isinstance(self._marker_set, pybropt.popgen.MarkerSet):
            raise RuntimeError("MarkerSet not provided")

    def calc_D(self):
        """
        Calculate linkage disequilibrium matrix as conventional D metric.

        Parameters
        ----------

        Returns
        -------
        D : numpy.ndarray
            A linkage disequilibrium matrix
        """
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
        sgeno = self._geno[:,sel,:]

        # calculate GEBVs
        gebv = self._genomic_model.predict(sgeno)

        # take the dot product if necessary
        if objcoeff is not None:
            gebv = gebv.dot(objcoeff)

        return gebv

    def lexsort(self, keys = None):
        """
        Get sorting indices using the marker set.
        """
        self._check_self_has_MarkerSet()

        # get indices
        indices = self._marker_set.lexsort(keys)

        # return indices
        return indices

    def reorder(self, indices):
        """
        Reorder the marker set and the genotype data.
        """
        # reorder genotypes
        self._geno = self._geno[:,:,indices]

        # reorder genetic model
        if self._genomic_model is not None:
            self._genomic_model.reorder(indices)

        # reorder marker set
        if self._marker_set is not None:
            self._marker_set.reorder(indices)

    def group(self):
        """
        Calculate grouping indices in the marker_set
        """
        # group marker set
        self._marker_set.group()

    # TODO: sort along taxa axis too
    def sort(self, keys = None):
        # get indices for sort
        indices = self.lexsort(keys)

        # reorder internals
        self.reorder(indices)

        # indicate that we've sorted
        self._sorted = True
        self._marker_set.sorted = True

        # calculate grouping indices
        self.group()

        if isinstance(self._marker_set, pybropt.popgen.GeneticMap):
            # remove any discrepancies in map
            self._marker_set.remove_discrepancies()

            # calculate crossover probabilities
            self._marker_set.calc_xo_prob()

    def remove(self, indices, axis, auto_sort = False):
        # if we have nothing, do nothing
        if indices is None:
            return

        # remove taxa names if we have them
        if (axis == 1) and (self._taxa is not None):
            self._taxa = numpy.delete(self._taxa, indices)

        # remove if we have a marker set
        if (axis == 2) and (self._marker_set is not None):
            self._marker_set.remove(indices, False) # do not auto sort

        # delete the keys
        self._geno = numpy.delete(self._geno, indices, axis = axis)

        # auto sort if needed
        if auto_sort:
            self.sort()

    def interpolate_genetic_map(self, base_genetic_map, kind = None, fill_value = None):
        """
        Interpolate the internal MarkerSet to a GeneticMap.
        """
        # check data types
        pybropt.util.check_is_GeneticMap(base_genetic_map, "base_genetic_map")

        # interpolate genetic map
        genetic_map = base_genetic_map.interpolate(
            self._marker_set.chr_grp,
            self._marker_set.chr_start,
            self._marker_set.chr_stop,
            self._marker_set.mkr_name,
            map_fncode = None,
            auto_sort = False,
            auto_mkr_rename = False,
            kind = kind,
            fill_value = fill_value
        )

        # replace marker set with genetic map
        self._marker_set = genetic_map

    def mkr_rename(self, new_mkr_name = None):
        self._check_self_has_MarkerSet()

        self._marker_set.mkr_rename(new_mkr_name)

    def mkr_mask(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None, invert = False):
        self._check_self_has_MarkerSet()

        # get mask
        mask = self._marker_set.mkr_mask(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            invert = invert
        )

        return mask

    def mkr_index(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None, invert = False):
        # find markers
        mask = self.mkr_mask(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            invert = invert
        )

        # get indices
        ix = numpy.flatnonzero(mask) if mask is not None else None

        return ix

    def mkr_remove(self, chr_grp = None, chr_start = None, chr_stop = None, mkr_name = None, invert = False, auto_sort = False):
        """
        Remove markers matching a unique signature from the marker set.
        """
        # get markers
        ix = self.mkr_index(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            invert = False
        )

        self.remove(ix, 2, auto_sort)

    def taxa_rename(self, new_taxa = None):
        if new_taxa is None:
            return

        self._taxa = new_taxa

    def taxa_mask(self, taxa = None, invert = False):
        # NOTE: I set this up as using tuples in case we need to grow
        masks = [None]

        # test whether self._taxa is in taxa
        if taxa is not None:
            masks[0] = numpy.in1d(self._taxa, taxa)

        # filter out None
        masks = list(m for m in masks if m is not None)

        # default value of None (for cases where len(masks) == 0)
        mask = None

        # if the len(masks) > 0, logical_and merge them together
        if len(masks) > 0:
            mask = numpy.logical_and.reduce(masks)

            # invert mask if we need to
            if invert:
                mask = ~mask

        # return mask
        return mask

    def taxa_index(self, taxa = None, invert = False):
        # find taxa
        mask = self.taxa_mask(
            taxa = taxa,
            invert = invert
        )

        # convert to indices
        ix = numpy.flatnonzero(mask) if mask is not None else None

        return ix

    # TODO: functionality for auto_sort of taxa names
    def taxa_remove(self, taxa = None, invert = False, auto_sort = False):
        # find taxa
        ix = self.taxa_index(
            taxa = taxa,
            invert = invert
        )

        # remove taxa; no auto_sort yet
        self.remove(ix, 1, False)

    ############################################################################
    ############################# Static Methods ###############################
    ############################################################################
    @staticmethod
    def from_array(geno, chr_grp, chr_start, chr_stop, mkr_name = None, taxa = None,
        base_genetic_map = None, genomic_model = None, auto_sort = False,
        auto_mkr_rename = False, kind = None, fill_value = None):
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
        base_genetic_map : GeneticMap
        taxa : numpy.ndarray

        Returns
        -------
        population : Population
            A Population object.
        """
        # check data types
        pybropt.util.check_is_matrix(geno, "geno")
        pybropt.util.check_is_matrix(chr_grp, "chr_grp")
        pybropt.util.check_is_matrix(chr_start, "chr_start")
        pybropt.util.check_is_matrix(chr_stop, "chr_stop")
        pybropt.util.cond_check_is_matrix(mkr_name, "mkr_name")
        pybropt.util.cond_check_is_GeneticMap(base_genetic_map, "base_genetic_map")
        pybropt.util.cond_check_is_GenomicModel(genomic_model, "genomic_model")
        pybropt.util.cond_check_matrix_dtype_is_string_(taxa, "taxa")

        # check dimensions
        pybropt.util.check_matrix_ndim(geno, "geno", 3)
        pybropt.util.check_matrix_ndim(chr_grp, "chr_grp", 1)
        pybropt.util.check_matrix_ndim(chr_start, "chr_start", 1)
        pybropt.util.check_matrix_ndim(chr_stop, "chr_stop", 1)
        pybropt.util.cond_check_matrix_ndim(mkr_name, "mkr_name", 1)
        pybropt.util.cond_check_matrix_ndim(taxa, "taxa", 1)

        # check matrix compatiblity lengths
        nloci = geno.shape[2]
        pybropt.util.check_matrix_size(chr_grp, "chr_grp", nloci)
        pybropt.util.check_matrix_size(chr_start, "chr_start", nloci)
        pybropt.util.check_matrix_size(chr_stop, "chr_stop", nloci)
        pybropt.util.cond_check_matrix_size(mkr_name, "mkr_name", nloci)
        pybropt.util.cond_check_matrix_size(taxa, "taxa", geno.shape[1])


        # make marker set; do not sort or rename
        marker_set = pybropt.popgen.MarkerSet(
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            auto_sort = False,
            auto_mkr_rename = False
        )

        # declare output variable
        population = Population(
            geno = geno,
            marker_set = marker_set,
            genomic_model = genomic_model,
            taxa = taxa
        )

        # interpolate genetic map positions using base_genetic_map
        population.interpolate_genetic_map(
            base_genetic_map = base_genetic_map,
            kind = kind,
            fill_value = fill_value
        )

        # sort if needed
        if auto_sort:
            population.sort()

        # rename markers if needed
        if auto_mkr_rename:
            population.marker_set.mkr_rename()

        return population

    @staticmethod
    def from_vcf(fname, base_genetic_map = None, genomic_model = None,
        auto_sort = False, auto_mkr_rename = False, kind = None, fill_value = None):
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
            pybropt.util.check_matrix_all_value(phases[:,2], "is_phased", True)

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

        population = Population.from_array(
            geno = geno,
            chr_grp = chr_grp,
            chr_start = chr_start,
            chr_stop = chr_stop,
            mkr_name = mkr_name,
            taxa = taxa,
            base_genetic_map = base_genetic_map,
            genomic_model = genomic_model,
            auto_sort = auto_sort,
            auto_mkr_rename = auto_mkr_rename,
            kind = kind,
            fill_value = fill_value
        )

        return population