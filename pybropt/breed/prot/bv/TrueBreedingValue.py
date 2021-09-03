from . import BreedingValueProtocol

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
    def estimate(self, obj, gmat, gpmod = None, **kwargs):
        """
        Estimate breeding values.

        Parameters
        ----------
        pt_or_bv : PhenotypeDataFrame, BreedingValueMatrix
            Phenotype dataframe or breeding value matrix to use to estimate
            breeding values.
            Unused by this class.
        gmat : PhasedGenotypeMatrix
            Genotype matrix to use for breeding value calculation.
            Also used to align genotypes in estimation output.

        Returns
        -------
        bvmat : BreedingValueMatrix
            Breeding value matrix.
        """
        # check inputs
        check_is_PhasedGenotypeMatrix(gmat, "gmat")

        # get default parameters
        if gpmod is None:
            gpmod = self.gpmod
        else:
            check_is_GenomicModel(gpmod, "gpmod")

        # calculate true breeding values
        bvmat = gpmod.predict(gmat)

        return bvmat
