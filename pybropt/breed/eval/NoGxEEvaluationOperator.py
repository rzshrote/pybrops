import numpy

from . import EvaluationOperator

from pybropt.popgen.bvmat import DenseEstimatedBreedingValueMatrix

class NoGxEEvaluationOperator(EvaluationOperator):
    """docstring for NoGxEEvaluationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nenv, var_E, rng, **kwargs):
        """
        nenv : int
            Number of environments.
        var_E : float or array_like of floats
            Environmental variance (constant) for each environment.
        rng : numpy.random.Generator
            Random number generator for introducing stochastic environmental
            error.
        """
        super(NoGxEEvaluationOperator, self).__init__(**kwargs)
        self.nenv = nenv
        self.var_E = var_E
        self.rng = rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def evaluate(self, t_cur, t_max, pgvmat, gmod, **kwargs):
        """
        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        pgvmat : PhasedGenotypeVariantMatrix
            Genotypes to evaluate.
        gmod_true : GenomicModel
            True genomic model.

        Returns
        -------
        out : tuple
            A tuple containing two elements: (bvmat, misc)
            bvmat : BreedingValueMatrix
                A matrix of breeding values
            bvmat_true : BreedingValueMatrix
                A matrix of true breeding values
            misc : dict
                Miscellaneous output (user defined).
        """
        # get true breeding values
        bvmat_true = gmod_true.pred(pgvmat)

        # generate raw phenotypes: (n,t) + (r,n,t) -> (r,n,t)
        raw = bvmat_true.mat + self.rng.normal(
            0.0,                                                # no deviation/GxE
            self.var_E,                                         # float or (t,)
            (self.nenv, bvmat_true.ntaxa, bvmat_true.ntrait),   # (r, n, t)
        )

        # calculate mean across all environments: (t,n,t) -> (n,t)
        mat = raw.mean(0)

        # create DenseEstimatedBreedingValueMatrix
        bvmat = DenseEstimatedBreedingValueMatrix(
            mat = mat,
            raw = raw,
            trait = bvmat_true.trait,
            taxa = bvmat_true.taxa,
            taxa_grp = bvmat_true.taxa_grp,
        )

        # preserve grouping metadata if it exists
        bvmat.taxa_grp_name = bvmat_true.taxa_grp_name
        bvmat.taxa_grp_stix = bvmat_true.taxa_grp_stix
        bvmat.taxa_grp_spix = bvmat_true.taxa_grp_spix
        bvmat.taxa_grp_len = bvmat_true.taxa_grp_len

        # empty dictionary of misc output
        misc = {}

        return bvmat, bvmat_true, misc



################################################################################
################################## Utilities ###################################
################################################################################
def is_NoGxEEvaluationOperator(v):
    return isinstance(v, NoGxEEvaluationOperator)

def check_is_NoGxEEvaluationOperator(v, vname):
    if not isinstance(v, NoGxEEvaluationOperator):
        raise TypeError("variable '{0}' must be a NoGxEEvaluationOperator".format(vname))

def cond_check_is_NoGxEEvaluationOperator(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_NoGxEEvaluationOperator(v, vname)
