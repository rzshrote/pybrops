from pybropt.core.mat import TaxaVariantMatrix

class HaplotypeMatrix(TaxaVariantMatrix):
    """
    An abstract class for haplotype matrix objects.

    The purpose of this abstract class is to proved base functionality for:
        1) haplotype metadata manipulation routines.
        2) haplotype math functions.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        HaplotypeMatrix constructor

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(HaplotypeMatrix, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## General matrix properties ###############
    def ploidy():
        doc = "Ploidy number represented by matrix property."
        def fget(self):
            """Get matrix ploidy number"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set matrix ploidy number"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete matrix ploidy number"""
            raise NotImplementedError("method is abstract")
        return locals()
    ploidy = property(**ploidy())

    # REMARK: this property is defined in PhasedMatrix as well. the purpose of
    # adding this here as well is to facilitate determination of the number of
    # phases represented by the HaplotypeMatrix. If 0, unphased; if >0, phased.
    def nphase():
        doc = "Number of chromosome phases represented by the matrix."
        def fget(self):
            """Get number of phases"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set number of phases"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete number of phases"""
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

    ############## Matrix summary statistics ###############
    def thcount(self, dtype):
        """
        Haplotype count of the non-zero haplotype within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, h) containing haplotype counts of the
            haplotype coded as 1 for all 'n' individuals, for all 'h'
            haplotypes.
        """
        raise NotImplementedError("method is abstract")

    def thfreq(self, dtype):
        """
        Haplotype frequency of the non-zero haplotype within each taxon.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array and of the accumulator in which the
            elements are summed.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (n, h) containing haplotype frequencies of
            the haplotype coded as 1 for all 'n' individuals, for all 'h'
            haplotypes.
        """
        raise NotImplementedError("method is abstract")

    def hcount(self, dtype):
        """
        Haplotype count of the non-zero haplotype across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype counts of the
            haplotype coded as 1 for all 'h' haplotypes.
        """
        raise NotImplementedError("method is abstract")

    def hfreq(self, dtype):
        """
        Haplotype frequency of the non-zero haplotype across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype frequencies of
            the haplotype coded as 1 for all 'h' haplotypes.
        """
        raise NotImplementedError("method is abstract")

    def mhf(self, dtype):
        """
        Minor haplotype frequency across all taxa.

        Parameters
        ----------
        dtype : dtype, optional
            The type of the returned array. If None, return native dtype.

        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape (h,) containing haplotype frequencies for
            the minor haplotype.
        """
        raise NotImplementedError("method is abstract")

    def mehe(self):
        """
        Mean expected heterozygosity across all taxa.

        Returns
        -------
        mehe : numpy.float64
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
            An int64 array of shape (3, h) containing haplotype counts across
            all 'h' haplotypes.
            Rows are as follows:
                out[0] = count of '0' genotype across all loci
                out[1] = count of '1' genotype across all loci
                out[2] = count of '2' genotype across all loci
        """
        raise NotImplementedError("method is abstract")

    def gtfreq(self):
        """
        Gather haplotype frequencies for homozygous major, heterozygous,
        homozygous minor across all individuals.

        Returns
        -------
        out : numpy.ndarray
            An float64 array of shape (3, h) containing haplotype counts across
            all 'h' haplotypes.
            Rows are as follows:
                out[0] = frequency of '0' genotype across all loci
                out[1] = frequency of '1' genotype across all loci
                out[2] = frequency of '2' genotype across all loci
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
def is_HaplotypeMatrix(v):
    return isinstance(v, HaplotypeMatrix)

def check_is_HaplotypeMatrix(v, varname):
    if not isinstance(v, HaplotypeMatrix):
        raise TypeError("'%s' must be a HaplotypeMatrix." % varname)

def cond_check_is_HaplotypeMatrix(v, varname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_HaplotypeMatrix(v, varname)