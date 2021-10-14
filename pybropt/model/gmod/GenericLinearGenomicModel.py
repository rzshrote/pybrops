import copy
import numpy
import h5py

from . import LinearGenomicModel

from pybropt.core.error import check_file_exists
from pybropt.core.error import check_group_in_hdf5
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
from pybropt.core.error import error_readonly

from pybropt.core.util import save_dict_to_hdf5
from pybropt.core.util import is_ndarray
from pybropt.popgen.gmat import is_GenotypeMatrix
from pybropt.popgen.bvmat import is_BreedingValueMatrix
from pybropt.popgen.bvmat import DenseGenomicEstimatedBreedingValueMatrix

class GenericLinearGenomicModel(LinearGenomicModel):
    """docstring for GenericLinearGenomicModel."""

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, beta, u, trait = None, model_name = None, params = None, **kwargs):
        """
        Constructor for GenericLinearGenomicModel class.

        Parameters
        ----------
        beta : numpy.ndarray
            A numpy.float64 fixed effect regression coefficient matrix of shape (q,t).
            Where:
                q : is the number of fixed effect predictors (e.g. environments)
                t : is the number of individuals
        u : numpy.ndarray
            A numpy.float64 random effect regression coefficient matrix of shape (p,t).
            Where:
                p : is the number of random effect predictors (e.g. genomic markers)
                t : is the number of individuals
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
        self.beta = beta
        self.u = u
        self.trait = trait
        self.model_name = model_name
        self.params = params

    def __copy__(self):
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : GenomicModel
        """
        out = self.__class__(
            beta = copy.copy(self.beta),
            u = copy.copy(self.u),
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
            beta = copy.deepcopy(self.beta),
            u = copy.deepcopy(self.u),
            trait = copy.deepcopy(self.trait),
            model_name = copy.deepcopy(self.model_name),
            params = copy.deepcopy(self.params)
        )

        return out

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ############## Linear Genomic Model Data ###############
    def beta():
        doc = "Fixed effect regression coefficients"
        def fget(self):
            """Get fixed effect regression coefficients"""
            return self._beta
        def fset(self, value):
            """Set fixed effect regression coefficients"""
            check_is_ndarray(value, "beta")
            check_ndarray_ndim(value, "beta", 2)
            check_ndarray_dtype_is_float64(value, "beta")
            self._beta = value
        def fdel(self):
            """Delete fixed effect regression coefficients"""
            del self._beta
        return locals()
    beta = property(**beta())

    def u():
        doc = "Random effect regression coefficients"
        def fget(self):
            """Get random effect regression coefficients"""
            return self._u
        def fset(self, value):
            """Set random effect regression coefficients"""
            check_is_ndarray(value, "u")
            check_ndarray_ndim(value, "u", 2)
            check_ndarray_dtype_is_float64(value, "u")
            self._u = value
        def fdel(self):
            """Delete random effect regression coefficients"""
            del self._u
        return locals()
    u = property(**u())

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
            self._trait = value
        def fdel(self):
            del self._trait
        return locals()
    trait = property(**trait())

    def ntrait():
        doc = "Number of traits predicted by the model"
        def fget(self):
            """Get the number of traits predicted by the model"""
            return len(self._trait)
        def fset(self, value):
            """Set the number of traits predicted by the model"""
            error_readonly("ntrait")
        def fdel(self):
            """Delete the number of traits predicted by the model"""
            error_readonly("ntrait")
        return locals()
    ntrait = property(**ntrait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ####### methods for model fitting and prediction #######
    def fit_numpy(self, Y, X, Z, **kwargs):
        """
        Fit the model.

        Parameters
        ----------
        Y : numpy.ndarray
            A phenotype matrix of shape (n,t).
        X : numpy.ndarray
            A covariate matrix of shape (n,q).
        Z : numpy.ndarray
            A genotypes matrix of shape (n,p).
        trait : numpy.ndarray
            A trait name array of shape (t,).
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise AttributeError("GenericLinearGenomicModel is read-only")

    def fit(self, ptobj, cvobj, gtobj, **kwargs):
        """
        Fit the model.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, PhenotypeDataFrame, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        trait : numpy.ndarray, None
            A trait name array of shape (t,).
        **kwargs : **dict
            Additional keyword arguments.
        """
        raise AttributeError("GenericLinearGenomicModel is read-only")

    ######## methods for estimated breeding values #########
    def predict_numpy(self, X, Z, **kwargs):
        """
        Predict breeding values.

        Remark: The difference between 'predict_numpy' and 'gebv_numpy' is that
        'predict_numpy' can incorporate other factors (e.g., fixed effects) to
        provide prediction estimates.

        Parameters
        ----------
        X : numpy.ndarray
            A matrix of covariates.
        Z : numpy.ndarray
            A matrix of genotype values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        Y_hat : numpy.ndarray
            A matrix of predicted breeding values.
        """
        # Y = Xβ + Zu
        Y_hat = (X @ self.beta) + (Z @ self.u)

        return Y_hat

    def predict(self, cvobj, gtobj, **kwargs):
        """
        Predict breeding values.

        Remark: The difference between 'predict' and 'gebv' is that 'predict'
        can incorporate other factors (e.g., fixed effects) to provide
        prediction estimates.

        Parameters
        ----------
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix,
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            Estimated breeding values.
        """
        # process cvobj
        if is_ndarray(cvobj):
            X = cvobj
        else:
            raise TypeError("accepted types are numpy.ndarray")

        # process gtobj
        if is_GenotypeMatrix(gtobj):
            Z = gtobj.mat_asformat("{0,1,2}")
            taxa = gtobj.taxa
            taxa_grp = gtobj.taxa_grp
        elif is_ndarray(gtobj):
            Z = gtobj
            taxa = None
            taxa_grp = None
        else:
            raise TypeError("accepted types are GenotypeMatrix, numpy.ndarray")

        # make predictions
        Y_hat = self.predict_numpy(X, Z, **kwargs)

        # create output breeding value matrix
        out = DenseGenomicEstimatedBreedingValueMatrix(
            mat = Y_hat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self.trait
        )

        return out

    def score_numpy(self, Y, X, Z, **kwargs):
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        Y : numpy.ndarray
            A matrix of phenotypes.
        X : numpy.ndarray
            A matrix of covariates.
        Z : numpy.ndarray
            A matrix of genotypes.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape (t,).
            Where:
                t : is the number of traits.
        """
        # TODO: array shape checks

        # calculate predictions
        # (n,q) @ (q,t) -> (n,t)
        # (n,p) @ (p,t) -> (n,t)
        # (n,t) + (n,t) -> (n,t)
        Y_hat = (X @ self.beta) + (Z @ self.u)

        # calculate sum of squares error
        # (n,t) - (n,t) -> (n,t)
        # (n,t)**2 -> (n,t)
        # (n,t).sum(0) -> (t,)
        SSE = ((Y - Y_hat)**2).sum(0)

        # calculate means for each trait
        # (n,t).mean(0) -> (t,)
        Y_mean = Y.mean(0)

        # calculate sum of squares total
        # (n,t) - (t,) -> (n,t) - (1,t)
        # (n,t) - (1,t) -> (n,t)
        # (n,t).sum(0) -> (t,)
        SST = ((Y - Y_mean)**2).sum(0)

        # calculate R**2
        # (t,) / (t,) -> (t,)
        # scalar - (t,) -> (t,)
        Rsq = (1.0 - SSE/SST)

        return Rsq

    def score(self, ptobj, cvobj, gtobj, **kwargs):
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, PhenotypeDataFrame, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape (t,).
            Where:
                t : is the number of traits.
        """
        # process ptobj
        if is_BreedingValueMatrix(ptobj):
            Y = ptobj.mat
        elif is_PhenotypeDataFrame(ptobj):
            raise RuntimeError("not implmented yet")
        elif is_ndarray(ptobj):
            Y = ptobj
        else:
            raise TypeError("must be BreedingValueMatrix, PhenotypeDataFrame, numpy.ndarray")

        # process cvobj
        if is_ndarray(cvobj):
            X = cvobj
        else:
            raise TypeError("must be numpy.ndarray")

        # process gtobj
        if is_GenotypeMatrix(gtobj):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif is_ndarray(gtobj):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, numpy.ndarray")

        # calculate coefficient of determination
        Rsq = self.score_numpy(Y, X, Z, **kwargs)

        return Rsq

    ######## methods for estimated breeding values #########
    def gebv_numpy(self, Z, **kwargs):
        """
        Calculate genomic estimated breeding values.

        Remark: The difference between 'predict_numpy' and 'gebv_numpy' is that
        'predict_numpy' can incorporate other factors (e.g., fixed effects) to
        provide prediction estimates.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotype values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        gebv_hat : numpy.ndarray
            A matrix of genomic estimated breeding values.
        """
        # Y = Zu
        gebv_hat = (Z @ self.u)

        return gebv_hat

    def gebv(self, gtobj, **kwargs):
        """
        Calculate genomic estimated breeding values.

        Remark: The difference between 'predict' and 'gebv' is that 'predict'
        can incorporate other factors (e.g., fixed effects) to provide
        prediction estimates.

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            Genomic estimated breeding values.
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            Z = gtobj.mat_asformat("{0,1,2}")
            taxa = gtobj.taxa
            taxa_grp = gtobj.taxa_grp
        elif is_ndarray(gtobj):
            Z = gtobj
            taxa = None
            taxa_grp = None
        else:
            raise TypeError("accepted types are GenotypeMatrix, numpy.ndarray")

        # make predictions
        gebv_hat = self.gebv_numpy(Z, **kwargs)

        # create output breeding value matrix
        out = DenseGenomicEstimatedBreedingValueMatrix(
            mat = gebv_hat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self.trait
        )

        return out

    ###### methods for population variance prediction ######
    def var_G_numpy(self, Z, **kwargs):
        """
        Calculate the population genetic variance.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypes.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # estimate breeding values (n,t)
        gebv = self.gebv_numpy(Z, **kwargs)

        # calculate variance
        # (n,t).var(0) -> (t,)
        out = gebv.var(0)

        return out

    def var_G(self, gtobj, **kwargs):
        """
        Calculate the population genetic variance.

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif is_ndarray(gtobj):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        out =self.var_G_numpy(Z, **kwargs)

        return out

    def var_A_numpy(self, Z, **kwargs):
        """
        Calculate the population additive genetic variance

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypes.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # estimate breeding values (n,t)
        gebv = self.gebv_numpy(Z, **kwargs)

        # calculate variance
        # (n,t).var(0) -> (t,)
        out = gebv.var(0)

        return out

    def var_A(self, gtobj, **kwargs):
        """
        Calculate the population additive genetic variance

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif is_ndarray(gtobj):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        out =self.var_A_numpy(Z, **kwargs)

        return out

    def var_a_numpy(self, p, ploidy, **kwargs):
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : int
            Ploidy of the species.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # change shape to (p,1)
        p = p[:,None]

        # calculate additive genic variance
        # (p,t)**2 * (p,1) * (p,1) -> (p,t)
        # (p,t).sum[0] -> (t,)
        # scalar * (t,) -> (t,)
        out = (ploidy**2.0) * ((self.u**2) * p * (1.0 - p)).sum(0)

        return out

    def var_a(self, gtobj, ploidy = None, **kwargs):
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        ploidy : int
            Ploidy of the species.
            If ploidy is None:
                If gtobj is a GenotypeMatrix:
                    Get ploidy from GenotypeMatrix.
                If gtobj is a numpy.ndarray:
                    Assumed to be 2 (diploid).
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif is_ndarray(gtobj):
            if ploidy is None:
                ploidy = 2
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate genic variance
        out = self.var_a_numpy(p, ploidy, **kwargs)

        return out

    def bulmer_numpy(self, Z, p, ploidy, **kwargs):
        """
        Calculate the Bulmer effect.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypes.
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : int
            Ploidy of the species.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        sigma_A = self.var_A_numpy(Z)           # calculate additive genetic variance
        sigma_a = self.var_a_numpy(p, ploidy)   # calculate additive genetic variance
        mask = (sigma_a == 0.0)                 # determine where division by zero occurs
        denom = sigma_a.copy()                  # copy array
        denom[mask] = 1.0                       # substitute non-zero value
        out = sigma_A / denom                   # calculate Bulmer effect
        out[mask] = numpy.nan                   # add NaN's (avoids div by zero warning)
        return out

    def bulmer(self, gtobj, ploidy = None, **kwargs):
        """
        Calculate the Bulmer effect.

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        ploidy : int
            Ploidy of the species.
            If ploidy is None:
                If gtobj is a GenotypeMatrix:
                    Get ploidy from GenotypeMatrix.
                If gtobj is a numpy.ndarray:
                    Assumed to be 2 (diploid).
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            Array of Bulmer effects for each trait. In the event that additive
            genic variance is zero, NaN's are produced.
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            Z = gtobj.mat_asformat("{0,1,2}")
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif is_ndarray(gtobj):
            Z = gtobj       # get genotypes
            if ploidy is None:                  # if ploidy not provided
                ploidy = 2                      # assume diploid
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate Bulmer effect
        out = self.bulmer_numpy(Z, p, ploidy, **kwargs)

        return out

    ############# methods for selection limits #############
    def usl_numpy(self, p, ploidy, **kwargs):
        """
        Calculate the upper selection limit for a population.

        Parameters
        ----------
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : int
            Ploidy of the species.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # reshape allele frequencies
        # (p,) -> (p,1)
        p = p[:,None]

        # get maximum attainable genotype
        # (p,t) ? (p,1) : (p,1) -> (p,t)
        uslgeno = numpy.where(
            self.u > 0.0,       # if the allele effect is positive
            p > 0.0,            # +allele: 1 if we have at least one +allele
            p >= 1.0            # -allele: 1 if we have fixation for -allele
        )

        # calculate usl value
        # scalar * (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = (float(ploidy) * self.u * uslgeno).sum(0)

        return out

    def usl(self, gtobj, ploidy = None, **kwargs):
        """
        Calculate the upper selection limit for a population.

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif is_ndarray(gtobj):
            if ploidy is None:
                ploidy = 2
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate genic variance
        out = self.usl_numpy(p, ploidy, **kwargs)

        return out

    def lsl_numpy(self, p, ploidy, **kwargs):
        """
        Calculate the lower selection limit for a population.

        Parameters
        ----------
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : int
            Ploidy of the species.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # reshape allele frequencies
        # (p,) -> (p,1)
        p = p[:,None]

        # get minimum attainable genotype
        # (p,t) ? (p,1) : (p,1) -> (p,t)
        lslgeno = numpy.where(
            self.u > 0.0,       # if the allele effect is positive
            p >= 1.0,           # +allele: 1 if we have fixation for +allele
            p > 0.0             # -allele: 1 if we have at least one -allele
        )

        # calculate lsl value
        # scalar * (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = (float(ploidy) * self.u * lslgeno).sum(0)

        return out

    def lsl(self, gtobj, ploidy = None, **kwargs):
        """
        Calculate the lower selection limit for a population.

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        **kwargs : **dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        # process gtobj
        if is_GenotypeMatrix(gtobj):
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif is_ndarray(gtobj):
            if ploidy is None:
                ploidy = 2
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate genic variance
        out = self.lsl_numpy(p, ploidy, **kwargs)

        return out

    ################### File I/O methods ###################
    @staticmethod
    def from_hdf5(filename, groupname = None):
        """
        Read GenotypeMatrix from an HDF5 file.

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
        check_file_exists(filename)                             # check file exists
        h5file = h5py.File(filename, "r")                       # open HDF5 in read only
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            check_group_in_hdf5(groupname, h5file, filename)    # check that group exists
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### check that we have all required fields
        required_fields = ["beta", "u"]                         # all required arguments
        for field in required_fields:                           # for each required field
            fieldname = groupname + field                       # concatenate base groupname and field
            check_group_in_hdf5(fieldname, h5file, filename)    # check that group exists
        ######################################################### read data
        data_dict = {                                           # output dictionary
            "beta": None,
            "u" : None,
            "trait": None,
            "model_name": None,
            "params": None
        }
        data_dict["beta"] = h5file[groupname + "beta"][()]      # read beta array
        data_dict["u"] = h5file[groupname + "u"][()]            # read u array
        fieldname = groupname + "trait"                         # construct "groupname/trait"
        if fieldname in h5file:                                 # if "groupname/trait" in hdf5
            data_dict["trait"] = h5file[fieldname][()]          # read trait array
            data_dict["trait"] = numpy.object_(                 # convert trait string from byte to utf-8
                [s.decode("utf-8") for s in data_dict["trait"]]
            )
        fieldname = groupname + "model_name"                    # construct "groupname/model_name"
        if fieldname in h5file:                                 # if "groupname/model_name" in hdf5
            data_dict["model_name"] = h5file[fieldname][()]     # read string (as bytes); convert to utf-8
            data_dict["model_name"] = data_dict["model_name"].decode("utf-8")
        fieldname = groupname + "params"                        # construct "groupname/params"
        if fieldname in h5file:                                 # if "groupname/params" in hdf5
            data_dict["params"] = {}                            # create empty dictionary
            view = h5file[fieldname]                            # get view of dataset
            for key in view.keys():                             # for each field
                data_dict["params"][key] = view[key][()]        # extract data
        ######################################################### read conclusion
        h5file.close()                                          # close file
        ######################################################### create object
        glgmod = GenericLinearGenomicModel(**data_dict)         # create object from read data
        return glgmod

    def to_hdf5(self, filename, groupname = None):
        """
        Write GenotypeMatrix to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If None, GenotypeMatrix is written to the base HDF5 group.
        """
        h5file = h5py.File(filename, "a")                       # open HDF5 in write mode
        ######################################################### process groupname argument
        if isinstance(groupname, str):                          # if we have a string
            if groupname[-1] != '/':                            # if last character in string is not '/'
                groupname += '/'                                # add '/' to end of string
        elif groupname is None:                                 # else if groupname is None
            groupname = ""                                      # empty string
        else:                                                   # else raise error
            raise TypeError("'groupname' must be of type str or None")
        ######################################################### populate HDF5 file
        data_dict = {                                           # data dictionary
            "beta": self.beta,
            "u": self.u,
            "trait": self.trait,
            "model_name": self.model_name,
            "params": self.params
        }
        save_dict_to_hdf5(h5file, groupname, data_dict)         # write data
        ######################################################### write conclusion
        h5file.close()                                          # close the file



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenericLinearGenomicModel(v):
    """
    Determine whether an object is a GenericLinearGenomicModel.

    Parameters
    ----------
    v : object
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GenericLinearGenomicModel object instance.
    """
    return isinstance(v, GenericLinearGenomicModel)

def check_is_GenericLinearGenomicModel(v, vname):
    """
    Check if object is of type GenericLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GenericLinearGenomicModel):
        raise TypeError("variable '{0}' must be a GenericLinearGenomicModel".format(vname))

def cond_check_is_GenericLinearGenomicModel(v, vname, cond=(lambda s: s is not None)):
    """
    Conditionally check if object is of type GenericLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    cond : function
        A function returning True/False for whether to test if is a GenericLinearGenomicModel.
    """
    if cond(v):
        check_is_GenericLinearGenomicModel(v, vname)
