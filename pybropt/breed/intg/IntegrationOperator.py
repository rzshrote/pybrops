class IntegrationOperator:
    """docstring for IntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        super(IntegrationOperator, self).__init__()

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def integrate(self, t_cur, t_max, pgvmat, bvmat, geno, bval):
        """
        Integrate genotype and phenotype data into geno and bval dictionaries.

        Parameters
        ----------
        pgvmat : PhasedGenotypeVariantMatrix
        bvmat : BreedingValueMatrix
        geno : dict
        bval : dict

        Returns
        -------
        out : tuple
            (geno_new, bval_new, misc)
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_IntegrationOperator(v):
    return isinstance(v, IntegrationOperator)

def check_is_IntegrationOperator(v, vname):
    if not isinstance(v, IntegrationOperator):
        raise TypeError("variable '{0}' must be a IntegrationOperator".format(vname))

def cond_check_is_IntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_IntegrationOperator(v, vname)