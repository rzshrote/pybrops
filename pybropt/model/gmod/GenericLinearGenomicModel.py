import copy
import numpy

from . import LinearGenomicModel

from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_axis_len
from pybropt.core.error import check_ndarray_dtype_is_float64
from pybropt.core.error import check_ndarray_ndim

from pybropt.core.error import cond_check_is_dict
from pybropt.core.error import cond_check_is_ndarray
from pybropt.core.error import cond_check_is_str
from pybropt.core.error import cond_check_ndarray_axis_len
from pybropt.core.error import cond_check_ndarray_dtype_is_object
from pybropt.core.error import cond_check_ndarray_ndim

from pybropt.popgen.bvmat import DenseGenotypicEstimatedBreedingValueMatrix

class GenericLinearGenomicModel(LinearGenomicModel):
    """docstring for GenericLinearGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, mu, beta, trait = None, model_name = None, params = None, **kwargs):
        """
        Constructor for GenericLinearGenomicModel class.

        Parameters
        ----------
        mu : numpy.ndarray
            A numpy.float64 array of shape (t, 1).
            Where:
                t : is the number of traits.
        beta : numpy.ndarray
            A numpy.float64 array of shape (p, t).
            Where:
                p : is the number of genomic loci.
                t : is the number of traits.
        trait : numpy.ndarray, None
            A numpy.object_ array of shape (t,).
            Where:
                t : is the number of traits.
        model_name : str, None
        params : dict, None
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GenericLinearGenomicModel, self).__init__(**kwargs)

        # set variables
        self.mu = mu
        self.beta = beta
        self.trait = trait
        self.model_name = model_name
        self.params = params

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Linear Genomic Model Data ###############
    def mu():
        doc = "The mu property."
        def fget(self):
            return self._mu
        def fset(self, value):
            check_is_ndarray(value, "mu")
            check_ndarray_ndim(value, "mu", 2)
            check_ndarray_dtype_is_float64(value, "mu")
            check_ndarray_axis_len(value, "mu", 1, 1) # shape = (t, 1)
            self._mu = value
        def fdel(self):
            del self._mu
        return locals()
    mu = property(**mu())

    def beta():
        doc = "The beta property."
        def fget(self):
            return self._beta
        def fset(self, value):
            check_is_ndarray(value, "beta")
            check_ndarray_ndim(value, "beta", 2)
            check_ndarray_dtype_is_float64(value, "beta")
            check_ndarray_axis_len(value, "beta", 1, self._mu.shape[0]) # shape = (p, t)
            self._beta = value
        def fdel(self):
            del self._beta
        return locals()
    beta = property(**beta())

    ################## Genomic Model Data ##################
    def model_name():
        doc = "The model_name property."
        def fget(self):
            return self._model_name
        def fset(self, value):
            cond_check_is_str(value, "model_name")
            self._model_name = value
        def fdel(self):
            del self._model_name
        return locals()
    model_name = property(**model_name())

    def params():
        doc = "The params property."
        def fget(self):
            return self._params
        def fset(self, value):
            cond_check_is_dict(value, "params")
            self._params = value
        def fdel(self):
            del self._params
        return locals()
    params = property(**params())

    def trait():
        doc = "The trait property."
        def fget(self):
            return self._trait
        def fset(self, value):
            cond_check_is_ndarray(value, "trait")
            cond_check_ndarray_ndim(value, "trait", 1)
            cond_check_ndarray_dtype_is_object(value, "trait")
            cond_check_ndarray_axis_len(value, "trait", 0, self._mu.shape[0])
            self._trait = value
        def fdel(self):
            del self._trait
        return locals()
    trait = property(**trait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    def __copy__(self):
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : GenomicModel
        """
        out = self.__class__(
            mu = copy.copy(self.mu),
            beta = copy.copy(self.beta),
            trait = copy.copy(self.trait),
            model_name = copy.copy(self.model_name),
            params = copy.copy(self.params)
        )

        return out

    def __deepcopy__(self, memo):
        """
        Make a deep copy of the GenomicModel.

        Parameters
        ----------
        memo : dict

        Returns
        -------
        out : GenomicModel
        """
        out = self.__class__(
            mu = copy.deepcopy(self.mu),
            beta = copy.deepcopy(self.beta),
            trait = copy.deepcopy(self.trait),
            model_name = copy.deepcopy(self.model_name),
            params = copy.deepcopy(self.params)
        )

        return out

    ################# methods for model fitting and prediction #################
    def fit(self, gmat, bvmat):
        """
        Fit the model

        Parameters
        ----------
        gmat : GenotypeMatrix
        bvmat : BreedingValueMatrix
        """
        raise RuntimeError("GenericLinearGenomicModel is read-only")

    def predict(self, gmat):
        """
        Predict breeding values.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        gebvmat : GenotypicEstimatedBreedingValueMatrix
        """
        # get genotypes as 0,1,2
        geno = gmat.tacount()

        # calculate GEBVs
        gebv = self.mu.T + (geno @ self.beta)

        # create DenseGenotypicEstimatedBreedingValueMatrix
        gebvmat = DenseGenotypicEstimatedBreedingValueMatrix(
            mat = gebv,
            raw = None,
            se = None,
            trait = self.trait,
            taxa = gmat.taxa,
            taxa_grp = gmat.taxa_grp
        )

        return gebvmat

    def score(self, gmat, bvmat):
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotypes from which to predict breeding values.
        bvmat : BreedingValueMatrix
            True breeding values from which to score prediction accuracy.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape (t,).
            Where:
                t : is the number of traits.
        """
        # calculate prediction matrix
        y_pred = self.mu.T + ((gmat.tacount()) @ self.beta)

        # get pointer to true breeding values
        y_true = bvmat.mat

        # calculate sum of squares error
        SSE = ((y_true - y_pred)**2).sum()

        # calculate sum of squares total
        SST = ((y_true - y_true.mean())**2).sum()

        # calculate R**2
        Rsq = (1.0 - SSE/SST)

        return Rsq

    ################ methods for population variance prediction ################
    def var_G(self, gmat):
        """
        Calculate the population genetic variance.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        out : numpy.ndarray
        """
        bvmat = self.predict(gmat)  # make genotype predictions
        out = bvmat.mat.var(0)      # get variance for each trait
        return out

    def var_A(self, gmat):
        """
        Calculate the population additive genetic variance

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        out : numpy.ndarray
        """
        # identical to var_G since this model is completely additive
        bvmat = self.predict(gmat)  # make genotype predictions
        out = bvmat.mat.var(0)      # get variance for each trait
        return out

    def var_a(self, gmat):
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        out : numpy.ndarray
        """
        p = gmat.afreq()[:,None]    # (p,1) get allele frequencies
        # (p,t)**2 * (p,1) * (p,1) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = 4.0 * ((self.beta**2) * p * (1.0 - p)).sum(0)
        return out

    def bulmer(self, gmat):
        """
        Calculate the Bulmer effect.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        out : numpy.ndarray
        """
        sigma_A = self.var_A(gmat)
        sigma_a = self.var_a(gmat)
        return sigma_A / sigma_a

    ####################### methods for selection limits #######################
    def usl(self, gmat):
        """
        Calculate the upper selection limit for a population.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        out : numpy.ndarray
        """
        p = gmat.afreq()[:,None]    # (p,1) get allele frequencies
        maxgeno = numpy.where(      # (p,t) get maximum attainable genotype
            self.beta > 0.0,
            p > 0.0,
            p == 1.0
        )

        ploidy = float(gmat.ploidy) # get ploidy

        # scalar * (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = (ploidy * self.beta * maxgeno).sum(0)

        # (t,) + (t,1).flatten -> (t,)
        out += self.mu.flatten()

        return out

    def lsl(self, gmat):
        """
        Calculate the lower selection limit for a population.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        out : numpy.ndarray
        """
        p = gmat.afreq()[:,None]    # (p,1) get allele frequencies
        mingeno = numpy.where(      # (p,t) get minimum attainable genotype
            self.beta > 0.0,
            p == 1.0,
            p > 0.0
        )

        ploidy = float(gmat.ploidy) # get ploidy

        # scalar * (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = (ploidy * self.beta * mingeno).sum(0)

        # (t,) + (t,1).flatten -> (t,)
        out += self.mu.flatten()

        return out



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenericLinearGenomicModel(v):
    return isinstance(v, GenericLinearGenomicModel)

def check_is_GenericLinearGenomicModel(v, vname):
    if not isinstance(v, GenericLinearGenomicModel):
        raise TypeError("variable '{0}' must be a GenericLinearGenomicModel".format(vname))

def cond_check_is_GenericLinearGenomicModel(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenericLinearGenomicModel(v, vname)
