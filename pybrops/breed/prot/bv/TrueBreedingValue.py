from pybrops.breed.prot.bv.BreedingValueProtocol import BreedingValueProtocol

class TrueBreedingValue(BreedingValueProtocol):
    """docstring for TrueBreedingValue."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gpmod, **kwargs):
        super(TrueBreedingValue, self).__init__(**kwargs)
        self.gpmod = gpmod

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############### Genomic Model Properties ###############
    def gpmod():
        doc = "Genomic prediction model."
        def fget(self):
            """Get genomic prediction model"""
            return self._gpmod
        def fset(self, value):
            """Set genomic prediction model"""
            check_is_GenomicModel(value, "gpmod")
            self._gpmod = value
        def fdel(self):
            """Delete genomic prediction model"""
            del self._gpmod
        return locals()
    gpmod = property(**gpmod())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def estimate(self, ptobj, gtobj, miscout = None, gpmod = None, **kwargs):
        """
        Estimate breeding values.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, PhenotypeDataFrame, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        miscout : dict, None
            Pointer to a dictionary for miscellaneous user defined output.
            If dict, write to dict (may overwrite previously defined fields).
            If None, user defined output is not calculated or stored.
        gpmod : GenomicModel
            Genomic model used for predicting genotypes.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            A matrix of breeding values.
        """
        # check inputs
        check_is_PhasedGenotypeMatrix(gmat, "gmat")

        # get default parameters
        if gpmod is None:
            gpmod = self.gpmod
        else:
            check_is_GenomicModel(gpmod, "gpmod")

        # calculate true breeding values
        bvmat = gpmod.gebv(gmat)

        return bvmat
