class GenomicModel:
    """docstring for GenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs):
        """
        Constructor for the GenomicModel class.

        Parameters
        ----------
        **kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GenomicModel, self).__init__(**kwargs)

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################## Genomic Model Data ##################
    def model_name():
        doc = "The model_name property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    model_name = property(**model_name())

    def params():
        doc = "The params property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    params = property(**params())

    def trait():
        doc = "The trait property."
        def fget(self):
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            raise NotImplementedError("method is abstract")
        def fdel(self):
            raise NotImplementedError("method is abstract")
        return locals()
    trait = property(**trait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ################# methods for model fitting and prediction #################
    def fit(self, gmat, bvmat):
        """
        Fit the model

        Parameters
        ----------
        gmat : GenotypeMatrix
        bvmat : BreedingValueMatrix
        """
        raise NotImplementedError("method is abstract")

    def predict(self, gmat):
        """
        Predict breeding values.

        Parameters
        ----------
        gmat : GenotypeMatrix

        Returns
        -------
        bvmat : BreedingValueMatrix
        """
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

################################################################################
################################## Utilities ###################################
################################################################################
def is_GenomicModel(v):
    return isinstance(v, GenomicModel)

def check_is_GenomicModel(v, vname):
    if not isinstance(v, GenomicModel):
        raise TypeError("variable '{0}' must be a GenomicModel".format(vname))

def cond_check_is_GenomicModel(v, vname, cond=(lambda s: s is not None)):
    if cond(v):
        check_is_GenomicModel(v, vname)
