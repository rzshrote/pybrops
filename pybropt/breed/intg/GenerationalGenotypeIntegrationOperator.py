from . import GenotypeIntegrationOperator

class GenerationalGenotypeIntegrationOperator(GenotypeIntegrationOperator):
    """docstring for GenerationalGenotypeIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gqlen, gwind, gmult, **kwargs):
        super(GenerationalGenotypeIntegrationOperator, self).__init__(**kwargs)
        self.gqlen = gqlen
        self.gwind = gwind
        self.gmult = gmult

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def gintegrate(self, t_cur, t_max, pgvmat, geno):
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
        # copy dictionaries
        geno_new = dict(geno)

        # duplicate queue lists to avoid pointer problems
        geno_new["queue"] = list(geno["queue"])

        # process genotype queue
        while len(geno_new["queue"]) < self.gqlen:  # while queue is too short
            geno_new["queue"].append(None)          # append None to length
        geno_new["queue"].append(pgvmat)            # add pgvmat to end of queue
        new_geno = geno_new["queue"].pop(0)         # pop new genotypes from queue
        # print("0:", new_geno.taxa_grp)
        # calculate the taxa_grp minimum threshold ('+ 1' is necessary!)
        taxa_min = (t_cur - (self.gqlen + self.gwind) + 1) * self.gmult

        # process genotype main
        # print("taxa_min:", taxa_min)
        # print("taxa_grp:", geno_new["main"].taxa_grp)
        mask = geno_new["main"].taxa_grp < taxa_min    # create genotype mask
        # print("1:",geno_new["main"].mat.shape)
        geno_new["main"].remove(mask, axis = 1)         # delete old taxa
        # print("2:",geno_new["main"].mat.shape)
        geno_new["main"].append(                        # add new taxa
            values = new_geno.mat,
            axis = 1,
            taxa = new_geno.taxa,
            taxa_grp = new_geno.taxa_grp,
        )
        # print("3:",geno_new["main"].mat.shape)

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