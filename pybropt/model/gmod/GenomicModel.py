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

    def __copy__(self):
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : GenomicModel
        """
        raise NotImplementedError("method is abstract")

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
        raise NotImplementedError("method is abstract")

    ####### methods for model fitting and prediction #######
    def fit_raw(self, Y, Z, **kwargs):
        """
        Fit the model.

        Parameters
        ----------
        Y : numpy.ndarray
            A matrix of phenotypes.
        Z : numpy.ndarray
            A matrix of genotypes.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def fit(self, bvmat, gmat, **kwargs):
        """
        Fit the model

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            A matrix of breeding values.
        gmat : GenotypeMatrix
            A matrix of genotype values.
        **kwargs
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def predict_raw(self, Z, **kwargs):
        """
        Predict breeding values.

        Remark: The difference between 'predict_raw' and 'gebv_raw' is that
        'predict_raw' can incorporate other factors (e.g., fixed effects) to
        provide prediction estimates.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotype values.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        Y_hat : numpy.ndarray
            A matrix of predicted breeding values.
        """
        raise NotImplementedError("method is abstract")

    def predict(self, gmat, **kwargs):
        """
        Predict breeding values.

        Remark: The difference between 'predict' and 'gebv' is that 'predict'
        can incorporate other factors (e.g., fixed effects) to provide
        prediction estimates.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotypes from which to predict breeding values.

        Returns
        -------
        bvmat_hat : BreedingValueMatrix
            Estimated breeding values.
        """
        raise NotImplementedError("method is abstract")

    def score_raw(self, Y, Z, **kwargs):
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        Y : numpy.ndarray
            A matrix of true phenotypes.
        Z : numpy.ndarray
            A matrix of genotypes.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape (t,).
            Where:
                t : is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    def score(self, bvmat, gmat, **kwargs):
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        bvmat : BreedingValueMatrix
            True breeding values from which to score prediction accuracy.
        gmat : GenotypeMatrix
            Genotypes from which to predict breeding values.
        **kwargs
            Additional keyword arguments.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape (t,).
            Where:
                t : is the number of traits.
        """
        raise NotImplementedError("method is abstract")

    ######## methods for estimated breeding values #########
    def gebv_raw(self, Z):
        """
        Calculate genomic estimated breeding values.

        Remark: The difference between 'predict_raw' and 'gebv_raw' is that
        'predict_raw' can incorporate other factors (e.g., fixed effects) to
        provide prediction estimates.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotype values.

        Returns
        -------
        gebv_hat : numpy.ndarray
            A matrix of genomic estimated breeding values.
        """
        raise NotImplementedError("method is abstract")

    def gebv(self, gmat):
        """
        Calculate genomic estimated breeding values.

        Remark: The difference between 'predict' and 'gebv' is that 'predict'
        can incorporate other factors (e.g., fixed effects) to provide
        prediction estimates.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotypes from which to predict genomic estimated breeding values.

        Returns
        -------
        gebvmat_hat : BreedingValueMatrix
            Genomic estimated breeding values.
        """
        raise NotImplementedError("method is abstract")

    ###### methods for population variance prediction ######
    def var_G_raw(self, Z):
        """
        Calculate the population genetic variance.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypes.

        Returns
        -------
        out : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")

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

    def var_A_raw(self, Z):
        """
        Calculate the population additive genetic variance

        Parameters
        ----------
        Z : numpy.ndarray

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

    def var_a_raw(self, Z):
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        Z : numpy.ndarray

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

    def bulmer_raw(self, Z):
        """
        Calculate the Bulmer effect.

        Parameters
        ----------
        Z : numpy.ndarray

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

    ############# methods for selection limits #############
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

    ################### File I/O methods ###################
    @staticmethod
    def from_hdf5(filename, groupname):
        """
        Read GenomicModel from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is read from base HDF5 group.

        Returns
        -------
        gmat : GenotypeMatrix
            A genotype matrix read from file.
        """
        raise NotImplementedError("method is abstract")

    def to_hdf5(self, filename, groupname):
        """
        Write GenomicModel to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.
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
