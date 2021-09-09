from . import DenseGenotypeMatrix
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim
from pybropt.core.error import check_ndarray_dtype
from pybropt.core.error import error_readonly

class DensePhasedGenotypeMatrix(DenseGenotypeMatrix):
    """docstring for DensePhasedGenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mat, **kwargs):
        """
        DensePhasedGenotypeMatrix object constructor.

        Parameters
        ----------
        mat : numpy.int8
            A int8 binary genotype matrix of shape (m, n, p).
            Encoding is: {0, 1}
            Where:
                'm' is the number of chromosome phases (2 for diploid, etc.).
                'n' is the number of individuals.
                'p' is the number of markers.
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(DensePhasedGenotypeMatrix, self).__init__(
            mat = mat,
            **kwargs
        )

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################### Data Properites ####################
    def mat():
        doc = "The genotype matrix."
        def fget(self):
            """
            Retrieve a pointer to the genotype matrix.

            Returns
            -------
            mat : numpy.ndarray
                Pointer to the genotype matrix.
            """
            return self._mat
        def fset(self, value):
            """
            Set the pointer to the genotype matrix. Perform error checks.

            Parameters
            ----------
            mat : numpy.ndarray
                A int8 binary genotype matrix of shape (m, n, p).
                Where:
                    'm' is the number of chromosome phases (2 for diploid, etc.).
                    'n' is the number of individuals.
                    'p' is the number of markers.
                Cannot be None.
            """
            # check input data types, then set variable
            check_is_ndarray(value, "mat")
            check_ndarray_dtype(value, "mat", 'int8')
            check_ndarray_ndim(value, "mat", 3)
            self._mat = value
        def fdel(self):
            """
            Remove the genotype matrix (self._mat) from scope.
            """
            del self._mat
        return locals()
    mat = property(**mat())

    ################# Metadata Properites ##################
    def ploidy():
        doc = "The ploidy level represented by the genotype matrix."
        def fget(self):
            return self._mat.shape[0]
        def fset(self, value):
            error_readonly("ploidy")
        def fdel(self):
            error_readonly("ploidy")
        return locals()
    ploidy = property(**ploidy())

    def nphase():
        doc = "The nphase property."
        def fget(self):
            return self._mat.shape[0]
        def fset(self, value):
            error_readonly("nphase")
        def fdel(self):
            error_readonly("nphase")
        return locals()
    nphase = property(**nphase())

    def ntaxa():
        doc = "The number of taxa represented by the genotype matrix."
        def fget(self):
            return self._mat.shape[1]
        def fset(self, value):
            error_readonly("ntaxa")
        def fdel(self):
            error_readonly("ntaxa")
        return locals()
    ntaxa = property(**ntaxa())

    def nloci():
        doc = "The number of marker loci represented by the genotype matrix."
        def fget(self):
            return self._mat.shape[2]
        def fset(self, value):
            error_readonly("nloci")
        def fdel(self):
            error_readonly("nloci")
        return locals()
    nloci = property(**nloci())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############## Matrix summary statistics ###############
    def tacount(self, dtype = None):
        """
        Allele count within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, p) containing allele counts of the
            allele coded as 1 for all 'n' individuals, for all 'p' loci.
        """
        out = self._mat.sum(0)              # take sum across the phase axis (0)
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def tafreq(self, dtype = None):
        """
        Allele frequency within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        tafreq : numpy.ndarray
            A numpy.ndarray of shape (n, p) containing allele frequencies of the
            allele coded as 1 for all 'n' individuals, for all 'p' loci.
        """
        rnphase = 1.0 / self._mat.shape[0]  # take the reciprocal of the number of phases = 1 / nphase
        out = rnphase * self._mat.sum(0)    # take sum across the phase axis (0) and divide by nphase
        if dtype is not None:               # if dtype is specified
            dtype = numpy.dtype(dtype)      # ensure conversion to dtype class
            if out.dtype != dtype:          # if output dtype and desired are different
                out = dtype.type(out)       # convert to correct dtype
        return out

    def acount(self):
        """
        Allele count across all taxa.

        Returns
        -------
        acount : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele counts of the allele
            coded as 1 for all 'p' loci.
        """
        # take sum across the phase (0) and taxa (1) axes.
        return self._mat.sum((0,1))

    def afreq(self):
        """
        Allele frequency across all taxa.

        Returns
        -------
        afreq : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele frequencies of the
            allele coded as 1 for all 'p' loci.
        """
        # take the reciprocal of the number of phases 1 / ploidy*ntaxa
        rnphase = 1.0 / (self._mat.shape[0] * self._mat.shape[1])
        # take sum across the phase (0) and taxa (1) axes; divide by ploidy*ntaxa
        return rnphase * self._mat.sum((0,1))

    def maf(self):
        """
        Minor allele frequency across all taxa.

        Returns
        -------
        maf : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele frequencies for the
            minor allele.
        """
        # get minor allele frequencies
        out = self.afreq()

        # create mask of allele frequencies greater than 0.5
        mask = out > 0.5

        # calculate 1.0 - allele frequency for frequencies greater than 0.5
        out[mask] = 1.0 - out[mask]

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
        p = self.afreq()            # get allele frequency (p)
        out = (p * (1.0 - p)).sum() # take p*(1-p) across all loci and sum the products
        s = self._mat.shape         # get matrix shape
        out *= (s[0] / s[2])        # multiply summation by (nphase / nloci)
        return numpy.float64(out)

    def gtcount(self):
        """
        Gather genotype counts across for homozygous major, heterozygous,
        homozygous minor all individuals.

        Returns
        -------
        out : numpy.ndarray
            An int64 array of shape (3, p) containing genotype counts across
            all loci.
            Rows are as follows:
                out[0] = count of '0' genotype across all loci
                out[1] = count of '1' genotype across all loci
                out[2] = count of '2' genotype across all loci
        """
        # get genotypes as {0, 1, 2}
        gt = self._mat.sum(0)

        # allocate output array
        out = numpy.empty((3,self._mat.shape[2]), dtype='int64')

        # get genotype counts
        out[0] = (gt == 0).sum(0)       # homozygous 0/0
        out[1] = (gt == 1).sum(0)       # heterozygous 0/1, 1/0
        out[2] = (gt == 2).sum(0)       # homozygous 1/1

        return out

    def gtfreq(self):
        """
        Gather genotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Returns
        -------
        out : numpy.ndarray
            An float64 array of shape (3, p) containing genotype counts across
            all loci.
            Rows are as follows:
                out[0] = frequency of '0' genotype across all loci
                out[1] = frequency of '1' genotype across all loci
                out[2] = frequency of '2' genotype across all loci
        """
        recip = 1.0 / self._mat.shape[1]    # get reciprocal of number of taxa
        out = recip * self.gtcount()        # calculate genotype frequencies
        return out


################################################################################
################################## Utilities ###################################
################################################################################
def is_DensePhasedGenotypeMatrix(v):
    return isinstance(v, DensePhasedGenotypeMatrix)

def check_is_DensePhasedGenotypeMatrix(v, varname):
    if not isinstance(v, DensePhasedGenotypeMatrix):
        raise TypeError("'%s' must be a DensePhasedGenotypeMatrix." % varname)

def cond_check_is_DensePhasedGenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_DensePhasedGenotypeMatrix(v, varname)
