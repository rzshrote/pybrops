from . import LinearGenomicModel
from pybropt.core.error import check_is_ndarray
from pybropt.core.error import check_ndarray_ndim

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
        beta : numpy.ndarray
        trait : numpy.ndarray, None
        model_name : str, None
        params : dict, None
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GenericLinearGenomicModel, self).__init__(
            mu = mu,
            beta = beta,
            trait = trait,
            model_name = model_name,
            params = params,
            **kwargs
        )
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
            pybropt.util.check_matrix_dtype(value, "mu", numpy.float64)
            pybropt.util.check_matrix_axis_len(value, "mu", 1, 1) # shape = (t, 1)
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
            pybropt.util.check_matrix_dtype(value, "beta", numpy.float64)
            pybropt.util.check_matrix_axis_len(value, "beta", 1, self._mu.shape[0]) # shape = (p, t)
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
            pybropt.util.cond_check_is_string(value, "model_name")
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
            pybropt.util.cond_check_is_dict(value, "params")
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
            pybropt.util.cond_check_is_matrix(value, "trait")
            pybropt.util.cond_check_matrix_ndim(value, "trait", 1)
            pybropt.util.cond_check_matrix_dtype(value, "trait", numpy.string_)
            pybropt.util.cond_check_matrix_axis_len(value, "trait", 0, self._mu.shape[0])
            self._trait = value
        def fdel(self):
            del self._trait
        return locals()
    trait = property(**trait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def pred(self, gmat):
        """
        Predict breeding values.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        gebvmat : GenotypicEstimatedBreedingValueMatrix
        """
        raise NotImplementedError("method needs to be implemented")

    def score(self, gmat, bvmat):
        """
        Score the model's accuracy.

        Parameters
        ----------
        gmat : GenotypeMatrix
        bvmat : BreedingValueMatrix
        """
        raise NotImplementedError("method needs to be implemented")
