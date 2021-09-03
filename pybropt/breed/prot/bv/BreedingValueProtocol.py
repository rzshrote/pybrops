class BreedingValueProtocol:
    """docstring for BreedingValueProtocol."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(BreedingValueProtocol, self).__init__()

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def estimate(self, obj, gmat, **kwargs):
        """
        Estimate breeding values.

        Parameters
        ----------
        obj : PhenotypeDataFrame, BreedingValueMatrix
            Phenotype dataframe or breeding value matrix to use to estimate
            breeding values.
        gmat : GenotypeMatrix
            Genotype matrix to use for estimation. Also used to align genotypes
            in estimation output.

        Returns
        -------
        bvmat : BreedingValueMatrix
            Breeding value matrix.
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_BreedingValueProtocol(v):
    return isinstance(v, BreedingValueProtocol)

def check_is_BreedingValueProtocol(v, vname):
    if not isinstance(v, BreedingValueProtocol):
        raise TypeError("variable '{0}' must be a BreedingValueProtocol".format(vname))

def cond_check_is_BreedingValueProtocol(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_BreedingValueProtocol(v, vname)
