from pybrops.core.mat.TaxaVariantMatrix import TaxaVariantMatrix
from pybrops.popgen.gmap.GeneticMappableMatrix import GeneticMappableMatrix
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput

class GenotypeMatrix(TaxaVariantMatrix,GeneticMappableMatrix,HDF5InputOutput):
    """docstring for GenotypeMatrix."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for abstract class GenotypeMatrix.

        Parameters
        ----------
        kwargs : dict
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
            """Get ploidy level"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set ploidy level"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete ploidy level"""
            raise NotImplementedError("method is abstract")
        return locals()
    ploidy = property(**ploidy())

    def nphase():
        doc = "The number of phases represented by the genotype matrix."
        def fget(self):
            """Get number of phases represented by the genotype matrix"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of phases represented by the genotype matrix"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of phases represented by the genotype matrix"""
            raise NotImplementedError("method is abstract")
        return locals()
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
        return locals()
    mat_format = property(**mat_format())

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
        raise NotImplementedError("method is abstract")

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
        Gather haplotype counts across for homozygous major, heterozygous,
        homozygous minor all individuals.

        Returns
        -------
        out : numpy.ndarray
            An ``int64`` array of shape ``(3,h)`` containing haplotype counts
            across all ``h`` haplotypes.

            Where:

            - ``out[0]`` is the count of ``0`` genotype across all loci
            - ``out[1]`` is the count of ``1`` genotype across all loci
            - ``out[2]`` is the count of ``2`` genotype across all loci
        """
        raise NotImplementedError("method is abstract")

    def gtfreq(self):
        """
        Gather haplotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Returns
        -------
        out : numpy.ndarray
            An ``float64`` array of shape ``(3,h)`` containing haplotype counts
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
def is_GenotypeMatrix(v):
    """Return whether an object is a GenotypeMatrix or not"""
    return isinstance(v, GenotypeMatrix)

def check_is_GenotypeMatrix(v, varname):
    """Raise TypeError if object is not a GenotypeMatrix"""
    if not isinstance(v, GenotypeMatrix):
        raise TypeError("'%s' must be a GenotypeMatrix." % varname)

def cond_check_is_GenotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    """If object is not None, raise TypeError if object is not a GenotypeMatrix"""
    if cond(v):
        check_is_GenotypeMatrix(v, varname)