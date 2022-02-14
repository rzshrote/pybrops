import numpy

from . import EvaluationOperator

from pybrops.core import random as pbo_rng
from pybrops.core.error import check_is_int
from pybrops.core.error import check_is_positive
from pybrops.core.error import check_is_iterable
from pybrops.core.error import cond_check_is_Generator
from pybrops.popgen.bvmat import DenseEstimatedBreedingValueMatrix

class NoGxEEvaluationOperator(EvaluationOperator):
    """docstring for NoGxEEvaluationOperator."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, nenv, var_E, rng = None, **kwargs):
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
        check_is_int(nenv, "nenv")
        # TODO: check var_E
        cond_check_is_Generator(rng, "rng")

        # assign variables
        self.nenv = nenv
        self.var_E = var_E
        self.rng = pbo_rng if rng is None else rng

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################ environment parameters ################
    def nenv():
        doc = "The nenv property."
        def fget(self):
            return self._nenv
        def fset(self, value):
            check_is_int(value, "nenv")         # type check
            check_is_positive(value, "nenv")    # value check
            self._nenv = value
        def fdel(self):
            del self._nenv
        return locals()
    nenv = property(**nenv())

    def var_E():
        doc = "The var_E property."
        def fget(self):
            return self._var_E
        def fset(self, value):
            if numpy.issubdtype(type(value), numpy.number):
                check_is_positive(value, "var_E")   # make sure is positive
                value = [value]                     # construct list of copies
            elif not (hasattr(value, "__iter__") and hasattr(value, "__len__")):
                raise ValueError("variable 'var_E' must be iterable and have a length")
            self._var_E = numpy.float64(value)      # set values
        def fdel(self):
            del self._var_E
        return locals()
    var_E = property(**var_E())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def set_h2(self, pgvmat, gmod_true, h2):
        """
        Set the narrow sense heritability for environments.

        Parameters
        ----------
        pgvmat : PhasedGenotypeMatrix
            Founder genotypes.
        gmod_true : GenomicModel
            True genomic model.
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        """
        x = pgvmat.tacount()    # (n,p)
        b = gmod_true.beta      # (p,t)
        y = x @ b               # (n,p) @ (p,t) -> (n,t)
        var_A = y.var(0)        # (n,t) -> (t,)
        # TODO: determine if we should use genetic or genic variance for var_E

        # calculate environmental variance
        # var_E = (1 - h2)/h2 * var_A - var_G
        # we assume var_G is zero, so var_E = (1 - h2)/h2 * var_A
        # scalar - (t,) -> (t,)
        # (t,) / (t,) -> (t,)
        # (t,) * (t,) -> (t,)
        self.var_E = (1.0 - h2) / h2 * var_A

    def set_H2(self, pgvmat, gmod_true, H2):
        """
        Set the broad sense heritability for environments.

        Parameters
        ----------
        pgvmat : PhasedGenotypeMatrix
            Founder genotypes.
        gmod_true : GenomicModel
            True genomic model.
        h2 : float, numpy.ndarray
            Narrow sense heritability.
        """
        raise NotImplementedError("method is abstract")

    def evaluate(self, t_cur, t_max, pgvmat, gmod_true, **kwargs):
        """
        Parameters
        ----------
        t_cur : int
            Current generation number.
        t_max : int
            Maximum (deadline) generation number.
        pgvmat : PhasedGenotypeMatrix
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

        # create phenotype matrix shape
        rawshape = (self.nenv, bvmat_true.ntaxa, bvmat_true.ntrait)

        # generate raw phenotypes: (n,t) + (r,n,t) -> (r,n,t)
        raw = bvmat_true.mat + self.rng.normal(
            0.0,                            # no deviation/GxE
            numpy.sqrt(self._var_E),        # (1,) or (t,)
            rawshape                        # (r, n, t)
        )

        # calculate mean across all environments: (r,n,t) -> (n,t)
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
    # def from_h2(gmat, lgmod, nenv, h2, rng = None, **kwargs):
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
    def from_h2(gmat, lgmod, nenv, h2, rng = None, **kwargs):
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
