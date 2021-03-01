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
    def var_E():
        doc = "The var_E property."
        def fget(self):
            return self._var_E
        def fset(self, value):
            # TODO: check value is positive
            self._var_E = value
            self._std_E = numpy.sqrt(self._var_E)
        def fdel(self):
            del self._var_E
            del self._std_E
        return locals()
    var_E = property(**var_E())

    def std_E():
        doc = "The std_E property."
        def fget(self):
            return self._std_E
        def fset(self, value):
            # TODO: check value is positive
            self._std_E = value
            self._var_E = numpy.square(self._std_E)
        def fdel(self):
            del self._std_E
            del self._var_E
        return locals()
    std_E = property(**std_E())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def evaluate(self, t_cur, t_max, pgvmat, gmod_true, **kwargs):
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
        bvmat_true = gmod_true.predict(pgvmat)

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

    ############################################################################
    ############################## Static Methods ##############################
    ############################################################################
    # TODO: H^2 calculations
    # @staticmethod
    # def from_h2(gmat, lgmod, nenv, h2, rng, **kwargs):
    #     """
    #     h2 : float or numpy.ndarray
    #         Narrow sense heritability of trait for single rep evaluation.
    #
    #     Returns
    #     -------
    #     out : NoGxEEvaluationOperator
    #     """
    #     # get allele frequencies
    #     # (p,) -> (p,1)
    #     p = gmat.afreq()[:,None]    # (p,1)
    #     q = 1.0 - p                 # (p,1)
    #
    #     # get marker effects
    #     alpha = lgmod.beta          # (p,t)
    #
    #     # calculate additive variance: sum(2*p*q*alpha)
    #     # (p,1)*(p,1)*(p,t) -> (p,t)
    #     # (p,t).sum[0] -> (t,)
    #     var_A = 2.0 * (p*q*(alpha**2)).sum(0)
    #
    #     # calculate environmental variance
    #     # var_E = (1 - h2)/h2 * var_A - var_G
    #     # we assume var_G is zero, so var_E = (1 - h2)/h2 * var_A
    #     # scalar - (t,) -> (t,)
    #     # (t,) / (t,) -> (t,)
    #     # (t,) * (t,) -> (t,)
    #     var_E = (1.0 - h2) / h2 * var_A
    #     # print(var_E)
    #     out = NoGxEEvaluationOperator(
    #         nenv = nenv,
    #         var_E = var_E,
    #         rng = rng,
    #         **kwargs
    #     )
    #
    #     return out
    @staticmethod
    def from_h2(gmat, lgmod, nenv, h2, rng, **kwargs):
        """
        h2 : float or numpy.ndarray
            Narrow sense heritability of trait for single rep evaluation.

        Returns
        -------
        out : NoGxEEvaluationOperator
        """
        x = gmat.tacount()      # (n,p)
        b = lgmod.beta          # (p,t)
        y = x @ b               # (n,t)
        var_A = y.var(0)        # (n,t) -> (t,)

        # calculate environmental variance
        # var_E = (1 - h2)/h2 * var_A - var_G
        # we assume var_G is zero, so var_E = (1 - h2)/h2 * var_A
        # scalar - (t,) -> (t,)
        # (t,) / (t,) -> (t,)
        # (t,) * (t,) -> (t,)
        var_E = (1.0 - h2) / h2 * var_A
        # print(var_E)
        out = NoGxEEvaluationOperator(
            nenv = nenv,
            var_E = var_E,
            rng = rng,
            **kwargs
        )

        return out


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
