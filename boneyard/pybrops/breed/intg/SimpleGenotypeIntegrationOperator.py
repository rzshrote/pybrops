from . import GenotypeIntegrationOperator

class SimpleGenotypeIntegrationOperator(GenotypeIntegrationOperator):
    """docstring for SimpleGenotypeIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict):
        super(SimpleGenotypeIntegrationOperator, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def gintegrate(self, t_cur, t_max, pgvmat, geno):
        """
        Integrate genotype into geno dictionary.

        Integration is simple.
            1) Add genotype matrix to end of queue.
            2) Pop first genotype matrix off the queue and put into geno["main"].

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

        Returns
        -------
        out : tuple
            (geno_new, misc)
        """
        # copy dictionaries
        geno_new = dict(geno)

        # duplicate queue lists to avoid pointer problems
        geno_new["queue"] = list(geno["queue"])

        # process genotype queue
        geno_new["queue"].append(pgvmat)                # add pgvmat to end of queue
        geno_new["main"] = geno_new["queue"].pop(0)     # pop new genotypes from queue

        # empty dictionary
        misc = {}

        return geno_new, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_SimpleGenotypeIntegrationOperator(v):
    return isinstance(v, SimpleGenotypeIntegrationOperator)

def check_is_SimpleGenotypeIntegrationOperator(v, vname):
    if not isinstance(v, SimpleGenotypeIntegrationOperator):
        raise TypeError("variable '{0}' must be a SimpleGenotypeIntegrationOperator".format(vname))

def cond_check_is_SimpleGenotypeIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_SimpleGenotypeIntegrationOperator(v, vname)
