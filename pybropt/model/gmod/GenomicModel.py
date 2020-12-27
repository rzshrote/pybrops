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
            return self._params
        def fset(self, value):
            self._params = value
        def fdel(self):
            del self._params
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
    def fit(self, gmat, bvmat):
        """
        Fit the model

        Parameters
        ----------
        gmat : GenotypeMatrix
        bvmat : BreedingValueMatrix
        """
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

    def score(self, gmat, bvmat):
        """
        Score the model's accuracy.

        Parameters
        ----------
        gmat : GenotypeMatrix
        bvmat : BreedingValueMatrix
        """
        raise NotImplementedError("method is abstract")
