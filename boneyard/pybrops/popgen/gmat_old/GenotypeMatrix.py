from pybrops.core.mat import MutableMatrix

class GenotypeMatrix(MutableMatrix):
    """docstring for GenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        GenotypeMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GenotypeMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Matrix Metadata Properites ##############
    def ploidy():
        doc = "The ploidy level represented by the genotype matrix."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    ploidy = property(**ploidy())

    def nphase():
        doc = "The number of phases represented by the genotype matrix."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    nphase = property(**nphase())

    def ntaxa():
        doc = "The number of taxa represented by the genotype matrix."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    ntaxa = property(**ntaxa())

    def nloci():
        doc = "The number of marker loci represented by the genotype matrix."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    nloci = property(**nloci())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ############## Matrix summary statistics ###############
    def tacount(self, dtype):
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
        raise NotImplementedError("method is abstract")

    def tafreq(self, dtype):
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
        raise NotImplementedError("method is abstract")

    def acount(self):
        """
        Allele count of the non-zero allele across all taxa.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele counts of the allele
            coded as 1 for all 'p' loci.
        """
        raise NotImplementedError("method is abstract")

    def afreq(self):
        """
        Allele frequency of the non-zero allele across all taxa.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele frequencies of the
            allele coded as 1 for all 'p' loci.
        """
        raise NotImplementedError("method is abstract")

    def maf(self):
        """
        Minor allele frequency across all taxa.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (p) containing allele frequencies for the
            minor allele.
        """
        raise NotImplementedError("method is abstract")

    def mehe(self):
        """
        Mean expected heterozygosity across all taxa.

        Returns
        -------
        out : numpy.float64
            A 64-bit floating point representing the mean expected
            heterozygosity.
        """
        raise NotImplementedError("method is abstract")

    def gtcount(self):
        """
        Genotype counts for homozygous major, heterozygous, homozygous minor.
        """
        raise NotImplementedError("method is abstract")

    def gtfreq(self):
        """
        Genotype frequency for homozygous major, heterozygous, homozygous minor.
        """
        raise NotImplementedError("method is abstract")

    ################### Matrix File I/O ####################
    @staticmethod
    def from_hdf5(filename, groupname):
        """
        Read GenotypeMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : GenotypeMatrix
            A genotype matrix read from file.
        """
        raise NotImplementedError("method is abstract")

    def to_hdf5(self, filename, groupname):
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.
        """
        raise NotImplementedError("method is abstract")


################################################################################
################################## Utilities ###################################
################################################################################
def is_GenotypeMatrix(v):
    return isinstance(v, GenotypeMatrix)

def check_is_GenotypeMatrix(v, varname):
    if not isinstance(v, GenotypeMatrix):
        raise TypeError("'%s' must be a GenotypeMatrix." % varname)

def cond_check_is_GenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenotypeMatrix(v, varname)