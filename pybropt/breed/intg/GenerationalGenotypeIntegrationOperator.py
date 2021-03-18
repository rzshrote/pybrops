from . import GenotypeIntegrationOperator

class GenerationalGenotypeIntegrationOperator(GenotypeIntegrationOperator):
    """docstring for GenerationalGenotypeIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gwind, **kwargs):
        super(GenerationalGenotypeIntegrationOperator, self).__init__(**kwargs)
        self.gwind = gwind

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def gintegrate(self, t_cur, t_max, pgvmat, geno, gwind = None):
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

        Returns
        -------
        out : tuple
            (geno_new, misc)
        """
        # get genotype window length
        if gwind is None:
            gwind = self.gwind

        # copy dictionaries
        geno_new = dict(geno)

        # duplicate queue lists to avoid pointer problems
        geno_new["queue"] = list(geno["queue"])

        # get pointer to queue for future use
        queue = geno_new["queue"]

        # append new genotype matrix to back of queue
        queue.append(pgvmat)

        # pop first genotype matrix on queue off and discard
        queue.pop(0)

        # concatenate genotype matrices together
        geno_new["main"] = queue[0].concat_taxa(queue[0:gwind])

        # empty dictionary
        misc = {}

        return geno_new, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenerationalGenotypeIntegrationOperator(v):
    return isinstance(v, GenerationalGenotypeIntegrationOperator)

def check_is_GenerationalGenotypeIntegrationOperator(v, vname):
    if not isinstance(v, GenerationalGenotypeIntegrationOperator):
        raise TypeError("variable '{0}' must be a GenerationalGenotypeIntegrationOperator".format(vname))

def cond_check_is_GenerationalGenotypeIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenerationalGenotypeIntegrationOperator(v, vname)
