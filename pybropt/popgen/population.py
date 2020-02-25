import numpy
import os, sys
import cyvcf2
sys.path.append(os.path.dirname(os.path.dirname(os.path.realpath(__file__))))
# import our libraries
from util.error_subroutines import *

class population:
    """
    Object representing a population.

    Contains genotype data

    Parameters
    ==========
    geno : numpy.ndarray, list-like, default = None
        A (M, N, L) genotype matrix where 'M' is the number of chromosome phases,
        'N' is the number of individuals, and 'L' is the number of markers.
        Restrictions:
            1)  This matrix should be a binary matrix. Non-binary values will
                cause error.
            2)  Multiplying this matrix by 'coeff' must result in an estimated
                trait value (e.g. GEBV).

    coeff : numpy.ndarray, list-like, default = None
        A (T, L) regression coefficients matrix where 'T' is the number of
        trait effect sets, and 'L' is the number of markers.
        Note:
            The attribute 'wcoeff' is calculated from coeff.
            Where:
                wcoeff = abs(coeff) / abs(coeff).sum() along each axis.

    taxa : numpy.ndarray, list-like, default = None
        A (N,) taxa name matrix where 'N' is the number of individuals. The taxa
        names must be strings (or convertable to strings).

    chr_grp : numpy.ndarray, list-like, default = None
        A (L,) chromosome group matrix where 'L' is the number of markers. The
        chromsome groups must be strings (or convertable to strings).

    chr_pos : numpy.ndarray, list-like, default = None
        A (2, L) marker position matrix where 'L' is the number of markers, the
        first row of the matrix is the chromosome start position of the marker
        (1 indexed, inclusive), and the second row of the matrix is the
        chromsome stop position of the marker (1 indexed, inclusive).

    trait : numpy.ndarray, list-like, default = None
        A (T,) trait name matrix where 'T' is the number of traits.
        Ideally, a trait name string should be formatted like so:
            <trait_name:est_method>
            Where:
                trait_name is the trait name.
                est_method is the regression estimation method.
                Example:
                    "yield:rrBLUP"

    pheno : numpy.ndarray, list-like, default = None
        A (T, N) phenotype matrix where 'T' is the nubmer of traits, and 'N' is
        the number of individuals (taxa).


    Attributes
    ==========
    geno : numpy.ndarray
        A (M, N, L) binary genotype matrix where 'M' is the number of chromosome
        phases, 'N' is the number of individuals, and 'L' is the number of
        markers.
        Notes:
            Internally, 'geno' is stored as a numpy.ndarray(dtype = 'int8')

    coeff : numpy.ndarray
        A (T, L) regression coefficients matrix where 'T' is the number of
        trait effect sets, and 'L' is the number of markers.
        Notes:
            Internally, 'coeff' is stored as a numpy.ndarray(dtype = 'float64')

    wcoeff : numpy.ndarray
        A (T, L) weight coefficient matrix where 'T' is the number of trait
        effect sets, and 'L' is the number of markers.
        Notes:
            Internally, 'wcoeff' is stored as a numpy.ndarray(dtype = 'float64')
            The attribute 'wcoeff' is calculated from coeff.
            Where:
                wcoeff = abs(coeff) / abs(coeff).sum() along each axis.

    taxa : numpy.ndarray, list-like, default = None
        A (N,) taxa name matrix where 'N' is the number of individuals. The taxa
        names must be strings (or convertable to strings).
        Notes:
            Internally, 'taxa' is stored as a numpy.ndarray(dtype = 'object')

    chr_grp : numpy.ndarray, list-like, default = None
        A (L,) chromosome group matrix where 'L' is the number of markers. The
        chromsome groups must be strings (or convertable to strings).
        Notes:
            Internally, 'chr_grp' is stored as a numpy.ndarray(dtype = 'object')

    chr_pos : numpy.ndarray, list-like, default = None
        A (2, L) marker position matrix where 'L' is the number of markers, the
        first row of the matrix is the chromosome start position of the marker
        (1 indexed, inclusive), and the second row of the matrix is the
        chromsome stop position of the marker (1 indexed, inclusive).
        Notes:
            Internally, 'chr_pos' is stored as a numpy.ndarray(dtype = 'int64')

    chr_start : numpy.ndarray
        A (1, L) marker start position matrix where 'L' is the number of
        markers. This is an alias for accessing the first row of 'chr_pos'.
        Notes:
            Internally stored as a numpy.ndarray(dtype = 'int64')

    chr_stop : numpy.ndarray
        A (1, L) marker stop position matrix where 'L' is the number of
        markers. This is an alias for accessing the second row of 'chr_pos'.
        Notes:
            Internally stored as a numpy.ndarray(dtype = 'int64')

    trait : numpy.ndarray, list-like, default = None
        A (T,) trait name matrix where 'T' is the number of traits.
        Ideally, a trait name string should be formatted like so:
            <trait_name:est_method>
            Where:
                trait_name is the trait name.
                est_method is the regression estimation method.
                Example:
                    "yield:rrBLUP"
        Notes:
            Internally, 'trait' is stored as a numpy.ndarray(dtype = 'object')

    pheno : numpy.ndarray, list-like, default = None
        A (T, N) phenotype matrix where 'T' is the nubmer of traits, and 'N' is
        the number of individuals (taxa).
        Notes:
            Internally, 'pheno' is stored as a numpy.ndarray(dtype = 'float64')

    nmarker : int
        Number of markers.

    nindiv : int
        Number of individuals in the population.

    nphase : int
        Number of chromsome phases (2 for diploids, etc.).

    ntrait : int
        Number of traits.


    Methods
    =======
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, geno = None, coeff = None, taxa = None, chr_grp = None,
                 chr_pos = None, trait = None, pheno = None):
        """
        Construct a population object.

        Parameters
        ==========
        geno : numpy.ndarray, list-like, default = None
            A (M, N, L) genotype matrix where 'M' is the number of chromosome
            phases, 'N' is the number of individuals, and 'L' is the number of
            markers.
            Restrictions:
                1)  This matrix should be a binary matrix. Non-binary values
                    will cause error.
                2)  Multiplying this matrix by 'coeff' must result in an
                    estimated trait value (e.g. GEBV).

        coeff : numpy.ndarray, list-like, default = None
            A (L, T) regression coefficients matrix where 'T' is the number of
            trait effect sets, and 'L' is the number of markers.
            Note:
                The attribute 'wcoeff' is calculated from coeff.
                Where:
                    wcoeff = abs(coeff) / abs(coeff).sum() along each axis.

        taxa : numpy.ndarray, list-like, default = None
            A (N,) taxa name matrix where 'N' is the number of individuals.
            The taxa names must be strings (or convertable to strings).

        chr_grp : numpy.ndarray, list-like, default = None
            A (L,) chromosome group matrix where 'L' is the number of markers.
            The chromsome groups must be strings (or convertable to strings).

        chr_pos : numpy.ndarray, list-like, default = None
            A (2, L) marker position matrix where 'L' is the number of markers,
            the first row of the matrix is the chromosome start position of the
            marker (1 indexed, inclusive), and the second row of the matrix is
            the chromsome stop position of the marker (1 indexed, inclusive).

        trait : numpy.ndarray, list-like, default = None
            A (T,) trait name matrix where 'T' is the number of traits.
            Ideally, a trait name string should be formatted like so:
                <trait_name:est_method>
                Where:
                    trait_name is the trait name.
                    est_method is the regression estimation method.
                    Example:
                        "yield:rrBLUP"

        pheno : numpy.ndarray, list-like, default = None
            A (T, N) phenotype matrix where 'T' is the nubmer of traits, and 'N'
            is the number of taxa.
        """
        ########################################################################
        # process 'geno'
        if geno is not None:
            geno = numpy.array(geno)        # convert to a numpy.ndarray
            check_matrix_dtype_is_integer(geno, 'geno') # check dtype
            check_matrix_ndim(geno, 'geno', 3)          # check ndim
            check_matrix_is_binary(geno, 'geno')        # check for binary
            geno = numpy.int8(geno)         # convert to int8 matrix
        self._geno = geno   # set a private variable _geno

        ########################################################################
        # process 'coeff'
        # make a temporary variable for wcoeff
        wcoeff = None
        if coeff is not None:
            # convert 'coeff' to a numpy.ndarray
            coeff = numpy.array(coeff)
            # make sure 'coeff' is a numeric type
            check_matrix_dtype_is_numeric(coeff, 'coeff')# check dtype
            check_matrix_ndim(coeff, 'coeff', 2)# check ndim
            # convert matrix to float64 matrix
            coeff = numpy.float64(coeff)
            # calculate wcoeff matrix
            wcoeff = numpy.absolute(coeff)  # take abs(coeff)
            wcoeff /= wcoeff.sum(0)         # divide each col by its col sum
        # set private variables _coeff and _wcoeff
        self._coeff = coeff
        self._wcoeff = wcoeff

        ########################################################################
        # process 'taxa'
        if taxa is not None:
            # cast as strings, put into object numpy.ndarray
            taxa = numpy.object_([str(e) for e in taxa])
        # set private variable _taxa to hold taxa names
        self._taxa = taxa

        ########################################################################
        # process chr_grp
        if chr_grp is not None:
            # cast as strings, put into object numpy.ndarray
            chr_grp = numpy.object_([str(e) for e in chr_grp])
        # set private variable _chr_grp to hold chr_grp
        self._chr_grp = chr_grp

        ########################################################################
        # process chr_pos
        if chr_pos is not None:
            # convert 'chr_pos' to a numpy.ndarray
            chr_pos = numpy.array(chr_pos)
            # make sure 'chr_pos' is an integer type
            check_matrix_dtype_is_integer(chr_pos, 'chr_pos')
            check_matrix_ndim(chr_pos, 'chr_pos', 2)
            check_matrix_axis_len(chr_pos, 'chr_pos', 0, 2)
            # convert chr_pos to int64 type.
            chr_pos = numpy.int64(chr_pos)
        # set private variable _chr_pos to hold positions
        self._chr_pos = chr_pos

        ########################################################################
        # process traits
        if trait is not None:
            # cast as strings, put into object numpy.ndarray
            trait = numpy.object_([str(e) for e in trait])
        # set private variable _trait to hold traits
        self._trait = trait

        ########################################################################
        # process phenotypes
        if pheno is not None:
            # convert pheno to numpy.ndarray
            pheno = numpy.ndarray(pheno)
            check_matrix_dtype_is_numeric(pheno, 'pheno')
            check_matrix_ndim(pheno, 'pheno', 2)
            # convert pheno to float64 type
            pheno = numpy.float64(pheno)
        # set private variable _pheno to hold phenotypes
        self._pheno = pheno

        ########################################################################
        # set internal boolean to indicate the internals need to be validated
        self._internals_valid = False

        # final sanity checks to make sure matrix sizes are compatible
        self.validate_internals()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    # property for genotype matrix
    def geno():
        doc = "The geno property. Handles get, set, del for geno matrix."
        def fget(self):
            return self._geno
        def fset(self, value):
            # convert 'value' to a numpy.ndarray
            value = numpy.array(value)
            # make sure geno is an integer matrix
            check_matrix_dtype_is_integer(value, 'geno')
            check_matrix_ndim(value, 'geno', 3)
            check_matrix_is_binary(value, 'geno')
            # convert value to int8 matrix
            value = numpy.int8(value)
            # assign private variables
            self._geno = value
            self._internals_valid = False
        def fdel(self):
            del self._geno      # delete
            self._geno = None   # and set to None
        return locals()
    geno = property(**geno())

    ############################################################################
    # property for coeff matrix
    def coeff():
        doc = "The coeff property. Handles get, set, del for coeff matrix."
        def fget(self):
            return self._coeff
        def fset(self, value):
            # convert 'coeff' to a numpy.ndarray
            value = numpy.array(value)
            # make sure 'coeff' is a numeric type
            check_matrix_dtype_is_numeric(value, 'coeff')
            # convert matrix to float64 matrix
            value = numpy.float64(value)
            # set private variable _coeff
            self._coeff = value
            self._internals_valid = False
        def fdel(self):
            del self._coeff
            self._coeff = None
        return locals()
    coeff = property(**coeff())

    ############################################################################
    # property for wcoeff matrix
    def wcoeff():
        doc = "The wcoeff property. Handles get, set, del for wcoeff matrix."
        def fget(self):
            return self._wcoeff
        def fset(self, value):
            # convert 'wcoeff' to a numpy.ndarray
            value = numpy.array(value)
            # make sure 'wcoeff' is a numeric type
            check_matrix_dtype_is_numeric(value, 'wcoeff')
            # convert matrix to float64 matrix
            value = numpy.float64(value)
            # set private variable _wcoeff
            self._wcoeff = value
            self._internals_valid = False
        def fdel(self):
            del self._wcoeff
            self._wcoeff = None
        return locals()
    wcoeff = property(**wcoeff())

    ############################################################################
    # property for taxa matrix
    def taxa():
        doc = "The taxa property. Handles get, set, del for taxa matrix."
        def fget(self):
            return self._taxa
        def fset(self, value):
            # cast as strings, put into object numpy.ndarray
            value = numpy.object_([str(e) for e in value])
            # set private variable _taxa
            self._taxa = value
            self._internals_valid = False
        def fdel(self):
            del self._taxa
            self._taxa = None
        return locals()
    taxa = property(**taxa())

    ############################################################################
    # property for chr_grp matrix
    def chr_grp():
        doc = "The chr_grp property. Handles get, set, del for chr_grp matrix."
        def fget(self):
            return self._chr_grp
        def fset(self, value):
            # cast as strings, put into object numpy.ndarray
            value = numpy.object_([str(e) for e in value])
            # set private variable _chr_grp
            self._chr_grp = value
            self._internals_valid = False
        def fdel(self):
            del self._chr_grp
        return locals()
    chr_grp = property(**chr_grp())

    ############################################################################
    # property for chr_pos matrix
    def chr_pos():
        doc = "The chr_pos property. Handles get, set, del for chr_pos matrix."
        def fget(self):
            return self._chr_pos
        def fset(self, value):
            # convert 'chr_pos' to a numpy.ndarray
            value = numpy.array(value)
            # make sure 'chr_pos' is an integer type
            check_matrix_dtype_is_integer(value, 'chr_pos')
            # convert chr_pos to int64 type.
            value = numpy.int64(value)
            # set private variable _chr_pos
            self._chr_pos = value
            self._internals_valid = False
        def fdel(self):
            del self._chr_pos
            self._chr_pos = None
        return locals()
    chr_pos = property(**chr_pos())

    ############################################################################
    # property for chr_start
    def chr_start():
        doc = "The chr_pos property. Handles get, set, del for chr_pos matrix."
        def fget(self):
            return self._chr_pos[0]
        def fset(self, value):
            # convert 'chr_start' to a numpy.ndarray
            value = numpy.array(value)
            # make sure 'chr_start' is an integer type
            check_matrix_dtype_is_integer(value, 'chr_start')
            # make sure dimension lengths are correct
            if value.shape != (self._chr_pos.shape[1],):
                raise ValueError("'chr_start' must have the shape %s" %
                    ( (self._nmarker,) )
                )
            # convert chr_pos to int64 type.
            value = numpy.int64(value)
            # set private variable _chr_pos
            self._chr_pos[0] = value
        def fdel(self):
            raise RuntimeError("Cannot delete chr_start.")
        return locals()
    chr_start = property(**chr_start())

    ############################################################################
    # property for chr_stop
    def chr_stop():
        doc = "The chr_pos property. Handles get, set, del for chr_pos matrix."
        def fget(self):
            return self._chr_pos[1]
        def fset(self, value):
            # convert 'chr_stop' to a numpy.ndarray
            value = numpy.array(value)
            # make sure 'chr_stop' is an integer type
            check_matrix_dtype_is_integer(value, 'chr_stop')
            # make sure dimension lengths are correct
            if value.shape != (self._nmarker,):
                raise ValueError("'chr_stop' must have the shape %s" %
                    ( (self._nmarker,) )
                )
            # convert chr_pos to int64 type.
            value = numpy.int64(value)
            # set private variable _chr_pos
            self._chr_pos[1] = value
        def fdel(self):
            raise RuntimeError("Cannot delete chr_stop.")
        return locals()
    chr_stop = property(**chr_stop())

    ############################################################################
    # property for trait matrix
    def trait():
        doc = "The trait property. Handles get, set, del for trait matrix."
        def fget(self):
            return self._trait
        def fset(self, value):
            # cast as strings, put into object numpy.ndarray
            value = numpy.object_([str(e) for e in value])
            # set private variable _trait
            self._trait = value
            self._internals_valid = False
        def fdel(self):
            del self._trait
            self._trait = None
        return locals()
    trait = property(**trait())

    ############################################################################
    # property for pheno matrix
    def pheno():
        doc = "The pheno property. Handles get, set, del for pheno matrix."
        def fget(self):
            return self._pheno
        def fset(self, value):
            # convert pheno to numpy.ndarray
            value = numpy.ndarray(value)
            # make sure pheno is a numeric type
            check_matrix_dtype_is_numeric(value, 'pheno')
            # convert pheno to float64 type
            value = numpy.float64(value)
            # set private variable _pheno
            self._pheno = value
            self._internals_valid = False
        def fdel(self):
            del self._pheno
            self._pheno = None
        return locals()
    pheno = property(**pheno())

    ############################################################################
    # property for nphase
    def nphase():
        doc = "The nphase property. Handles get, set, del for nphase."
        def fget(self):
            return self._nphase
        def fset(self, value):
            raise ValueError("'nphase' is read-only.")
        def fdel(self):
            raise ValueError("'nphase' is read-only.")
        return locals()
    nphase = property(**nphase())

    ############################################################################
    # property for nindiv
    def nindiv():
        doc = "The nindiv property. Handles get, set, del for nindiv."
        def fget(self):
            return self._nindiv
        def fset(self, value):
            raise ValueError("'nindiv' is read-only.")
        def fdel(self):
            raise ValueError("'nindiv' is read-only.")
        return locals()
    nindiv = property(**nindiv())

    ############################################################################
    # property for nmarker
    def nmarker():
        doc = "The nmarker property. Handles get, set, del for nmarker."
        def fget(self):
            return self._nmarker
        def fset(self, value):
            raise ValueError("'nmarker' is read-only.")
        def fdel(self):
            raise ValueError("'nmarker' is read-only.")
        return locals()
    nmarker = property(**nmarker())

    ############################################################################
    # property for ntrait
    def ntrait():
        doc = "The ntrait property. Handles get, set, del for ntaxa."
        def fget(self):
            return self._ntrait
        def fset(self, value):
            raise ValueError("'ntrait' is read-only.")
        def fdel(self):
            raise ValueError("'ntrait' is read-only.")
        return locals()
    ntrait = property(**ntrait())

    ############################################################################
    ########################## Object Static Methods ###########################
    ############################################################################
    # geno = None, coeff = None, taxa = None, chr_grp = None,
    # chr_pos = None, trait = None, pheno = None

    @staticmethod
    def from_vcf(fname):
        # make VCF iterator
        vcf = cyvcf2.VCF(fname)

        # extract taxa names from vcf header
        taxa = numpy.object_(vcf.samples)

        # make empty lists to store extracted values
        geno = []
        chr_grp = []
        chr_pos = []

        # iterate through VCF file and accumulate variants
        for variant in vcf:
            # append chromosome string
            chr_grp.append(variant.CHROM)

            # append variant position coordinates
            chr_pos.append([variant.POS, variant.POS + len(variant.REF) - 1])

            # extract allele states + whether they are phased or not
            phases = numpy.int8(variant.genotypes)

            # check that they are all phased
            check_matrix_all_value(phases[:,2], "is_phased", True)

            # append genotype states
            geno.append(phases[:,0:2].copy())

        # convert and transpose genotype matrix
        geno = numpy.int8(geno).transpose(2,1,0).copy()

        # convert to object array
        chr_grp = numpy.object_(chr_grp)

        # convert and transpose chr_pos matrix
        chr_pos = numpy.int64(chr_pos).T.copy()

        # make output object
        outpop = population(
            geno = geno,
            taxa = taxa,
            chr_grp = chr_grp,
            chr_pos = chr_pos
        )

        return outpop

    ############################################################################
    ########################## Object Class Methods ############################
    ############################################################################
    @classmethod
    def validate_internals(self):
        # if we are internally valid
        if self._internals_valid:
            return

        # begin testing matrix shapes
        if self._geno is not None:
            M, N, L = self._geno.shape
            # test L dimensions
            if self._coeff is not None:
                if L != self._coeff.shape[0]:
                    raise ValueError("'geno' and 'coeff' are incompatible.")
            if self._wcoeff is not None:
                if L != self._wcoeff.shape[0]:
                    raise ValueError("'geno' and 'wcoeff' are incompatible.")
            if self._chr_grp is not None:
                if L != self._chr_grp.shape[0]:
                    raise ValueError("'geno' and 'chr_grp' are incompatible.")
            if self._chr_pos is not None:
                if L != self._chr_pos.shape[1]:
                    raise ValueError("'geno' and 'chr_pos' are incompatible.")
            # test N dimensions
            if self._taxa is not None:
                if N != self._taxa.shape[0]:
                    raise ValueError("'geno' and 'taxa' are incompatible.")
            if self._pheno is not None:
                if N != self._pheno.shape[1]
                    raise ValueError("'geno' and 'pheno' are incompatible.")
        if self._coeff is not None:
            L, T = self._coeff.shape
            # test T dimensions
            if self._trait is not None:
                if T != self._trait.shape[0]:
                    raise ValueError("'coeff' and 'trait' are incompatible.")
            if self._pheno is not None:
                if T != self._pheno.shape[0]:
                    raise ValueError("'coeff' and 'trait' are incompatible.")
            # test L dimensions
            if self._chr_grp is not None:
                if L != self._chr_grp.shape[0]:
                    raise ValueError("'coeff' and 'chr_grp' are incompatible.")
            if self._chr_pos is not None:
                if L != self._chr_pos.shape[1]:
                    raise ValueError("'coeff' and 'chr_pos' are incompatible.")
            # test wcoeff
            if self._wcoeff is not None:
                if self._coeff.shape != self._wcoeff.shape:
                    raise ValueError("'coeff' and 'wcoeff' are incompatible.")
        if self._wcoeff is not None:
            L, T = self._wcoeff.shape
            # test T dimensions
            if self._trait is not None:
                if T != self._trait.shape[0]:
                    raise ValueError("'wcoeff' and 'trait' are incompatible.")
            if self._pheno is not None:
                if T != self._pheno.shape[0]:
                    raise ValueError("'wcoeff' and 'trait' are incompatible.")
            # test L dimensions
            if self._chr_grp is not None:
                if L != self._chr_grp.shape[0]:
                    raise ValueError("'wcoeff' and 'chr_grp' are incompatible.")
            if self._chr_pos is not None:
                if L != self._chr_pos.shape[1]:
                    raise ValueError("'wcoeff' and 'chr_pos' are incompatible.")
        if self._chr_grp is not None:
            L, = self._chr_grp.shape
            # test L dimensions
            if self._chr_pos is not None:
                if L != self._chr_pos.shape[1]:
                    raise ValueError("'chr_grp' and 'chr_pos' are incompatible.")
        if self._trait is not None:
            T, = self._trait.shape
            # test T dimensions
            if self._pheno is not None:
                if T != self._pheno.shape[1]:
                    raise ValueError("'trait' and 'pheno' are incompatible.")
        if self._taxa is not None:
            N, = self._taxa.shape
            # test N dimensions
            if self._pheno is not None:
                if N != self._pheno.shape[0]:
                    raise ValueError("'taxa' and 'pheno' are incompatible.")

        # attempt to get number of phases, individuals, and markers from geno.
        if self._geno is not None:
            self._nphase, self._nindiv, self._nmarker = geno.shape
        # attempt to get number of traits and number of markers from coeff
        if self._coeff is not None:
            self._nmarker, self._ntrait = self._coeff.shape
        # attempt to get number of traits and number of markers from wcoeff
        if self._wcoeff is not None:
            self._nmarker, self._ntrait = self._wcoeff.shape
        # attempt to get number taxa from taxa
        if self._taxa is not None:
            self._nindiv, = self._taxa.shape # need that comma before =
        # attempt to get number of markers from chr_grp
        if self._chr_grp is not None:
            self._nmarker, = self._chr_grp.shape # need that comma before =
        # attempt to get number of markers from chr_pos
        if self._chr_pos is not None:
            self._nmarker = self._chr_pos.shape[1] # only want 2nd element
        # attempt to get number of traits from trait
        if self._trait is not None:
            self._ntrait, = self._trait.shape # need that comma before =
        # attempt to get number of traits and number of taxa from pheno
        if self._pheno is not None:
            self._ntrait, self._nindiv = self._pheno.shape

        return
