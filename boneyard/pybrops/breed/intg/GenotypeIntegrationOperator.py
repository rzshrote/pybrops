class GenotypeIntegrationOperator:
    """docstring for GenotypeIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(GenotypeIntegrationOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def gintegrate(self, t_cur, t_max, pgvmat, geno, **kwargs: dict):
        """
        Integrate genotype into geno dictionary.

        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        pgvmat : PhasedGenotypeVariantMatrix
            Genotype matrix to integrate.
        geno : dict
            Genotype dictionary into which to integrate.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        out : tuple
            (geno_new, misc)
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenotypeIntegrationOperator(v):
    return isinstance(v, GenotypeIntegrationOperator)

def check_is_GenotypeIntegrationOperator(v, vname):
    if not isinstance(v, GenotypeIntegrationOperator):
        raise TypeError("variable '{0}' must be a GenotypeIntegrationOperator".format(vname))

def cond_check_is_GenotypeIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenotypeIntegrationOperator(v, vname)
