from . import IntegrationOperator

class GenerationalIntegrationOperator(IntegrationOperator):
    """docstring for GenerationalIntegrationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, gqlen, gwind, gmult, **kwargs):
        """
        Parameters
        ----------
        gqlen : int
            Genotype queue length.
        gwind : int
            Generation window. Number of generations to keep in main breeding
            pool.
        gmult : int
            Generation multiplier. Used to differentiate between generational
            cohorts.
        """
        super(GenerationalIntegrationOperator, self).__init__(**kwargs)
        self.gqlen = gqlen
        self.gwind = gwind
        self.gmult = gmult

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def integrate(self, t_cur, t_max, pgvmat, bvmat, bvmat_true, geno, bval):
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
        # copy dictionaries
        geno_new = dict(geno)
        bval_new = dict(bval)

        # duplicate queue lists to avoid pointer problems
        geno_new["queue"] = list(geno["queue"])
        bval_new["queue"] = list(bval["queue"])

        # process genotype queue
        while len(geno_new["queue"]) < self.gqlen:  # while queue is too short
            geno_new["queue"].append(None)          # append None to length
        geno_new["queue"].append(pgvmat)            # add pgvmat to end of queue
        new_geno = geno_new["queue"].pop(0)         # pop new genotypes from queue

        # process breeding value queue
        while len(bval_new["queue"]) < self.gqlen:  # while queue is too short
            bval_new["queue"].append(None)          # append None to length
        bval_new["queue"].append(bvmat)             # add bvmat to end of queue
        new_bval = bval_new["queue"].pop(0)         # pop new breeding values from queue

        # process breeding value queue_true
        while len(bval_new["queue_true"]) < self.gqlen:  # while queue is too short
            bval_new["queue_true"].append(None)          # append None to length
        bval_new["queue_true"].append(bvmat)             # add bvmat to end of queue
        new_bval_true = bval_new["queue_true"].pop(0)    # pop new breeding values from queue

        # calculate the taxa_grp minimum threshold
        taxa_min = (t_cur - (self.gqlen + self.gwind)) * gmult

        # process genotype main
        mask = geno_new["main"].taxa_grp < taxa_min     # create genotype mask
        geno_new["main"].delete(mask, axis = 1)         # delete old taxa
        geno_new["main"].append(                        # add new taxa
            values = new_geno.mat,
            axis = 1,
            taxa = new_geno.taxa,
            taxa_grp = new_geno.taxa_grp,
        )

        # process breeding value matrix
        mask = bval_new["main"].taxa_grp < taxa_min     # create breeding value mask
        bval_new["main"].delete(mask, axis = 1)         # delete old taxa
        bval_new["main"].append(                        # add new taxa
            values = new_bval.mat,
            axis = 1,
            taxa = new_bval.taxa,
            taxa_grp = new_bval.taxa_grp
        )

        # process breeding value matrix
        mask = bval_new["main_true"].taxa_grp < taxa_min     # create breeding value mask
        bval_new["main_true"].delete(mask, axis = 1)         # delete old taxa
        bval_new["main_true"].append(                        # add new taxa
            values = new_bval.mat,
            axis = 1,
            taxa = new_bval.taxa,
            taxa_grp = new_bval.taxa_grp
        )

        misc = {}

        return geno_new, bval_new, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenerationalIntegrationOperator(v):
    return isinstance(v, GenerationalIntegrationOperator)

def check_is_GenerationalIntegrationOperator(v, vname):
    if not isinstance(v, GenerationalIntegrationOperator):
        raise TypeError("variable '{0}' must be a GenerationalIntegrationOperator".format(vname))

def cond_check_is_GenerationalIntegrationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenerationalIntegrationOperator(v, vname)