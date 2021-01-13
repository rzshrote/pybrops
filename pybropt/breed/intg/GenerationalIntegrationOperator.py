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
        # process genotype queue
        while len(geno["queue"]) < self.gqlen:  # while queue is too short
            geno["queue"].append(None)          # append None to length
        geno["queue"].append(pgvmat)            # add pgvmat to end of queue
        new_geno = geno["queue"].pop(0)         # pop new genotypes from queue

        # process breeding value queue
        while len(bval["queue"]) < self.gqlen:  # while queue is too short
            bval["queue"].append(None)          # append None to length
        bval["queue"].append(bvmat)             # add bvmat to end of queue
        new_bval = bval["queue"].pop(0)         # pop new breeding values from queue

        # calculate the taxa_grp minimum threshold
        taxa_min = (t_cur - (self.gqlen + self.gwind)) * gmult

        # process genotype main
        mask = geno["main"].taxa_grp < taxa_min     # create genotype mask
        geno["main"].delete(mask, axis = 1)         # delete old taxa
        geno["main"].append(                        # add new taxa
            values = new_geno.mat,
            axis = 1,
            taxa = None,
            taxa_grp = None,
        )


        raise NotImplementedError("method is abstract")



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
