"""
Module implementing classes and error checking routines for genomic prediction
models that incorporate genomic additive effects.
"""

import copy
from numbers import Integral
from typing import Dict, Optional, Sequence, Union
import h5py
import numpy
import pandas
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_io_h5py import check_group_in_hdf5
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_float64
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_eq, check_ndarray_ndim
from pybrops.core.error.error_type_python import check_is_dict, check_is_str_or_Sequence
from pybrops.core.error.error_type_python import check_is_str
from pybrops.core.error.error_type_numpy import check_ndarray_dtype_is_object
from pybrops.core.error.error_attr_python import error_readonly
from pybrops.core.error.error_value_python import check_dict_has_keys, check_len, check_str_value
from pybrops.core.io.CSVDictInputOutput import CSVDictInputOutput
from pybrops.core.io.PandasDictInputOutput import PandasDictInputOutput
from pybrops.core.util.h5py import save_dict_to_hdf5
from pybrops.model.gmod.AdditiveLinearGenomicModel import AdditiveLinearGenomicModel
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.bvmat.DenseGenomicEstimatedBreedingValueMatrix import DenseGenomicEstimatedBreedingValueMatrix

class DenseAdditiveLinearGenomicModel(
        AdditiveLinearGenomicModel,
        PandasDictInputOutput,
        CSVDictInputOutput,
    ):
    """
    The DenseAdditiveLinearGenomicModel class represents a Multivariate Multiple
    Linear Regression model.

    A Multivariate Multiple Linear Regression model is defined as:

    .. math::
        \\mathbf{Y} = \\mathbf{XB} + \\mathbf{ZU} + \\mathbf{E}

    Where:

    - :math:`\\mathbf{Y}` is a matrix of response variables of shape ``(n,t)``.
    - :math:`\\mathbf{X}` is a matrix of fixed effect predictors of shape ``(n,q)``.
    - :math:`\\mathbf{B}` is a matrix of fixed effect regression coefficients of shape ``(q,t)``.
    - :math:`\\mathbf{Z}` is a matrix of random effect predictors of shape ``(n,p)``.
    - :math:`\\mathbf{U}` is a matrix of random effect regression coefficients of shape ``(p,t)``.
    - :math:`\\mathbf{E}` is a matrix of error terms of shape ``(n,t)``.

    Block matrix modifications to :

    :math:`\\mathbf{Z}` and :math:`\\mathbf{U}` can be decomposed into block
    matrices pertaining to different sets of effects:

    .. math::
        \\mathbf{Z} = \\begin{bmatrix} \\mathbf{Z_{misc}} & \\mathbf{Z_{a}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{Z_{misc}}` is a matrix of miscellaneous random effect predictors of shape ``(n,p_misc)``
    - :math:`\\mathbf{Z_{a}}` is a matrix of additive genomic marker predictors of shape ``(n,p_a)``

    .. math::
        \\mathbf{U} = \\begin{bmatrix} \\mathbf{U_{misc}} \\\\ \\mathbf{U_{a}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{U_{misc}}` is a matrix of miscellaneous random effects of shape ``(p_misc,t)``
    - :math:`\\mathbf{U_{a}}` is a matrix of additive genomic marker effects of shape ``(p_a,t)``

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors.
    - ``p_misc`` is the number of miscellaneous random effect predictors.
    - ``p_a`` is the number of additive genomic marker predictors.
    - The sum of ``p_misc`` and ``p_a`` equals ``p``.
    - ``t`` is the number of traits
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self, 
            beta: numpy.ndarray, 
            u_misc: Union[numpy.ndarray,None], 
            u_a: Union[numpy.ndarray,None], 
            trait: Optional[numpy.ndarray] = None, 
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseAdditiveLinearGenomicModel class.

        Parameters
        ----------
        beta : numpy.ndarray
            A ``float64`` fixed effect regression coefficient matrix of shape
            ``(q,t)``.

            Where:

            - ``q`` is the number of fixed effect predictors (e.g. environments).
            - ``t`` is the number of traits.
        
        u_misc : numpy.ndarray, None
            A ``float64`` random effect regression coefficient matrix of shape
            ``(p_misc,t)`` containing miscellaneous effects.

            Where:

            - ``p_misc`` is the number of miscellaneous random effect predictors.
            - ``t`` is the number of traits.

            If ``None``, then set to an empty array of shape ``(0,t)``.
        
        u_a : numpy.ndarray, None
            A ``float64`` random effect regression coefficient matrix of shape
            ``(p_a,t)`` containing additive marker effects.

            Where:

            - ``p_a`` is the number of additive marker effect predictors.
            - ``t`` is the number of traits.

            If ``None``, then set to an empty array of shape ``(0,t)``.
        
        trait : numpy.ndarray, None
            An ``object_`` array of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        
        model_name : str, None
            Name of the model.
        
        hyperparams : dict, None
            Model parameters.
        
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        # set variables (order dependent)
        self.beta = beta
        self.u_misc = u_misc
        self.u_a = u_a
        self.trait = trait
        self.model_name = model_name
        self.hyperparams = hyperparams

    def __repr__(
            self
        ) -> str:
        """
        Return repr(self).
        
        Returns
        -------
        out : str
            A representation of the object.
        """
        return "<{0} of shape [B = {1}, U_misc = {2}, U_a = {3}] at {4}>".format(
            type(self).__name__,
            self.beta.shape,
            self.u_misc.shape,
            self.u_a.shape,
            hex(id(self))
        )

    def __copy__(
            self
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : DenseAdditiveLinearGenomicModel
            A shallow copy of the model.
        """
        out = self.__class__(
            beta       = copy.copy(self.beta),
            u_misc     = copy.copy(self.u_misc),
            u_a        = copy.copy(self.u_a),
            trait      = copy.copy(self.trait),
            model_name = copy.copy(self.model_name),
            hyperparams     = copy.copy(self.hyperparams),
        )

        return out

    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Make a deep copy of the GenomicModel.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseAdditiveLinearGenomicModel
            A deep copy of the model.
        """
        out = self.__class__(
            beta       = copy.deepcopy(self.beta, memo),
            u_misc     = copy.deepcopy(self.u_misc, memo),
            u_a        = copy.deepcopy(self.u_a, memo),
            trait      = copy.deepcopy(self.trait, memo),
            model_name = copy.deepcopy(self.model_name, memo),
            hyperparams     = copy.deepcopy(self.hyperparams, memo),
        )

        return out

    ############################ Object Properties #############################

    ############### Genomic Model Parameters ###############
    @property
    def nexplan(self) -> Integral:
        """Number of explanatory variables required by the model."""
        return self.nexplan_beta + self.nexplan_u

    @property
    def nparam(self) -> Integral:
        """Number of model parameters."""
        return self.nparam_beta + self.nparam_u

    ########### Linear Genomic Model Parameters ############
    @property
    def nexplan_beta(self) -> Integral:
        """Number of fixed effect explanatory variables required by the model."""
        return self._beta.shape[0]

    @property
    def nparam_beta(self) -> Integral:
        """Number of fixed effect parameters."""
        return self._beta.size

    @property
    def beta(self) -> numpy.ndarray:
        """Fixed effect regression coefficients."""
        return self._beta
    @beta.setter
    def beta(self, value: numpy.ndarray) -> None:
        """Set fixed effect regression coefficients"""
        check_is_ndarray(value, "beta")
        check_ndarray_ndim(value, "beta", 2)
        check_ndarray_dtype_is_float64(value, "beta")
        self._beta = value

    @property
    def nexplan_u(self) -> Integral:
        """Number of random effect explanatory variables required by the model."""
        return self.nexplan_u_misc + self.nexplan_u_a

    @property
    def nparam_u(self) -> Integral:
        """Number of random effect parameters."""
        return self.nparam_u_misc + self.nparam_u_a

    @property
    def u(self) -> numpy.ndarray:
        """Random effect regression coefficients."""
        # get random effects and concatenate along compatible axes
        out = numpy.concatenate([self.u_misc, self.u_a], axis = 0)
        return out
    @u.setter
    def u(self, value: numpy.ndarray) -> None:
        """Set random effect regression coefficients"""
        raise AttributeError("variable 'u' is read-only; use 'u_misc' and 'u_a' to modify 'u'.")

    ####### Additive Linear Genomic Model Parameters #######
    @property
    def nexplan_u_misc(self) -> Integral:
        """Number of miscellaneous random effect explanatory variables required by the model."""
        return self._u_misc.shape[0]

    @property
    def nparam_u_misc(self) -> Integral:
        """Number of miscellaneous random effect parameters."""
        return self._u_misc.size

    @property
    def u_misc(self) -> numpy.ndarray:
        """Miscellaneous random effect regression coefficients."""
        return self._u_misc
    @u_misc.setter
    def u_misc(self, value: numpy.ndarray) -> None:
        """Set miscellaneous random effect regression coefficients"""
        if value is None:                                   # if value is None
            t = self.beta.shape[1]                          # get number of traits
            value = numpy.empty((0,t), dtype = "float64")   # make empty array of shape (0,t)
        check_is_ndarray(value, "u_misc")
        check_ndarray_ndim(value, "u_misc", 2)
        check_ndarray_dtype_is_float64(value, "u_misc")
        self._u_misc = value

    @property
    def nexplan_u_a(self) -> Integral:
        """Number of additive genomic marker explanatory variables required by the model."""
        return self._u_a.shape[0]

    @property
    def nparam_u_a(self) -> Integral:
        """Number of additive genomic marker parameters."""
        return self._u_a.size

    @property
    def u_a(self) -> numpy.ndarray:
        """Additive genomic marker effects."""
        return self._u_a
    @u_a.setter
    def u_a(self, value: numpy.ndarray) -> None:
        """Set additive genomic marker effect regression coefficients"""
        if value is None:                                   # if value is None
            t = self.beta.shape[1]                          # get number of traits
            value = numpy.empty((0,t), dtype = "float64")   # make empty array of shape (0,t)
        check_is_ndarray(value, "u_a")
        check_ndarray_ndim(value, "u_a", 2)
        check_ndarray_dtype_is_float64(value, "u_a")
        self._u_a = value

    ################## Genomic Model Data ##################
    @property
    def model_name(self) -> str:
        """Description for property model_name."""
        return self._model_name
    @model_name.setter
    def model_name(self, value: Union[str,None]) -> None:
        """Set data for property model_name."""
        if value is None:
            value = ""
        check_is_str(value, "model_name")
        self._model_name = value
    
    @property
    def hyperparams(self) -> dict:
        """Description for property hyperparams."""
        return self._params
    @hyperparams.setter
    def hyperparams(self, value: Union[dict,None]) -> None:
        """Set data for property hyperparams."""
        if value is None:
            value = {}
        check_is_dict(value, "hyperparams")
        self._params = value
    
    @property
    def trait(self) -> numpy.ndarray:
        """Description for property trait."""
        return self._trait
    @trait.setter
    def trait(self, value: numpy.ndarray) -> None:
        """Set data for property trait."""
        if value is not None:
            check_is_ndarray(value, "trait")
            check_ndarray_ndim(value, "trait", 1)
            check_ndarray_dtype_is_object(value, "trait")
        self._trait = value
    
    @property
    def ntrait(self) -> int:
        """Number of traits predicted by the model."""
        return self._beta.shape[1]
    @ntrait.setter
    def ntrait(self, value: int) -> None:
        """Set the number of traits predicted by the model"""
        error_readonly("ntrait")

    ############################## Object Methods ##############################

    #################### Model copying #####################
    def copy(
            self
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : GenomicModel
            A shallow copy of the original GenomicModel
        """
        return self.__copy__()

    def deepcopy(
            self,
            memo: Optional[dict] = None
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Make a deep copy of the GenomicModel.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : GenomicModel
            A deep copy of the original GenomicModel
        """
        return self.__deepcopy__(memo)

    ####### methods for model fitting and prediction #######
    def fit_numpy(
            self, 
            Y: numpy.ndarray, 
            X: numpy.ndarray, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Fit a dense, additive linear genomic model.

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
        kwargs : dict
            Additional keyword arguments.
        """
        raise AttributeError("DenseAdditiveLinearGenomicModel is read-only")

    def fit(
            self, 
            ptobj: object, 
            cvobj: numpy.ndarray, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
        """
        Fit a dense, additive linear genomic model.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, pandas.DataFrame, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        trait : numpy.ndarray, None
            A trait name array of shape (t,).
        kwargs : dict
            Additional keyword arguments.
        """
        raise AttributeError("DenseAdditiveLinearGenomicModel is read-only")

    ######## methods for estimated breeding values #########
    def predict_numpy(
            self, 
            X: numpy.ndarray, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        Y_hat : numpy.ndarray
            A matrix of estimated breeding values.
        """
        # Y = Xβ + Zu
        Y_hat = (X @ self.beta) + (Z @ self.u)

        return Y_hat

    def predict(
            self, 
            cvobj: numpy.ndarray, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> BreedingValueMatrix:
        """
        Predict breeding values.

        Remark: The difference between 'predict' and 'gebv' is that 'predict'
        can incorporate other factors (e.g., fixed effects) to provide
        prediction estimates.

        Parameters
        ----------
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            Estimated breeding values matrix.
        """
        # process cvobj
        if isinstance(cvobj, numpy.ndarray):
            X = cvobj
        else:
            raise TypeError("accepted types are numpy.ndarray")

        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
            taxa = gtobj.taxa
            taxa_grp = gtobj.taxa_grp
        elif isinstance(gtobj, numpy.ndarray):
            Z = gtobj
            taxa = None
            taxa_grp = None
        else:
            raise TypeError("accepted types are GenotypeMatrix, numpy.ndarray")

        # make predictions
        Y_hat = self.predict_numpy(X, Z, **kwargs)

        # create output breeding value matrix
        out = DenseGenomicEstimatedBreedingValueMatrix.from_numpy(
            mat = Y_hat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self.trait
        )

        return out

    def score_numpy(
            self, 
            Y: numpy.ndarray, 
            X: numpy.ndarray, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
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

    def score(
            self, 
            ptobj: Union[BreedingValueMatrix,numpy.ndarray], 
            cvobj: numpy.ndarray, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        ptobj : BreedingValueMatrix, pandas.DataFrame, numpy.ndarray
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        Rsq : numpy.ndarray
            A coefficient of determination array of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # process ptobj
        if isinstance(ptobj, BreedingValueMatrix):
            Y = ptobj.unscale()
        elif isinstance(ptobj, pandas.DataFrame):
            raise RuntimeError("not implmented yet")
        elif isinstance(ptobj, numpy.ndarray):
            Y = ptobj
        else:
            raise TypeError("must be BreedingValueMatrix, pandas.DataFrame, numpy.ndarray")

        # process cvobj
        if isinstance(cvobj, numpy.ndarray):
            X = cvobj
        else:
            raise TypeError("must be numpy.ndarray")

        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif isinstance(gtobj, numpy.ndarray):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, numpy.ndarray")

        # calculate coefficient of determination
        Rsq = self.score_numpy(Y, X, Z, **kwargs)

        return Rsq

    ######## methods for estimated breeding values #########
    def gebv_numpy(
            self, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate genomic estimated breeding values.

        Remark: The difference between 'predict_numpy' and 'gebv_numpy' is that
        'predict_numpy' can incorporate other factors (e.g., fixed effects) to
        provide prediction estimates.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotype values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        gebv_hat : numpy.ndarray
            A matrix of genomic estimated breeding values.
        """
        # type checks
        check_is_ndarray(Z, "Z")
        check_ndarray_ndim(Z, "Z", 2)
        check_ndarray_axis_len_eq(Z, "Z", 1, self.nexplan_u_a)

        # Y = Zu
        gebv_hat = (Z @ self.u_a)

        return gebv_hat

    def gebv(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> BreedingValueMatrix:
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : BreedingValueMatrix
            Genomic estimated breeding values matrix.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
            taxa = gtobj.taxa
            taxa_grp = gtobj.taxa_grp
        elif isinstance(gtobj, numpy.ndarray):
            Z = gtobj
            taxa = None
            taxa_grp = None
        else:
            raise TypeError("accepted types are GenotypeMatrix, numpy.ndarray")

        # make predictions
        gebv_hat = self.gebv_numpy(Z, **kwargs)

        # construct contrast (X*) to take cell means assuming that the first
        # column in X is a corner (first fixed effect is the intercept)
        # construct a contrast to calculate the GEBV center
        nfixed = self.beta.shape[0]
        Xstar = numpy.empty((1,nfixed), dtype = self.beta.dtype)
        Xstar[0,0] = 1
        Xstar[0,1:] = 1/nfixed

        # calculate intercepts (location)
        location = Xstar @ self.beta

        # add location to gebv_hat
        gebv_hat += location

        # create output breeding value matrix
        out = DenseGenomicEstimatedBreedingValueMatrix.from_numpy(
            mat = gebv_hat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self.trait
        )

        return out

    ######## methods for estimated genotypic value #########
    def gegv_numpy(
            self,
            Z: numpy.ndarray,
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate genomic estimated genotypic values.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypic markers.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of genomic estimated genotypic values.
        """
        return self.gebv_numpy(Z = Z, **kwargs)
    
    def gegv(
            self,
            gtobj: GenotypeMatrix,
            **kwargs: dict
        ) -> BreedingValueMatrix:
        """
        Calculate genomic estimated genotypic values.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypic markers.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of genomic estimated genotypic values.
        """
        return self.gebv(gtobj = gtobj, **kwargs)

    ###### methods for population variance prediction ######
    def var_G_numpy(
            self, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the population genetic variance.

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypes.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population genetic variances.

            Where:

            - ``t`` is the number of traits.
        """
        # estimate breeding values (n,t)
        gebv = self.gebv_numpy(Z, **kwargs)

        # calculate variance
        # (n,t).var(0) -> (t,)
        out = gebv.var(0)

        return out

    def var_G(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the population genetic variance.

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population genetic variances.

            Where:

            - ``t`` is the number of traits.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif isinstance(gtobj, numpy.ndarray):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        out = self.var_G_numpy(Z, **kwargs)

        return out

    def var_A_numpy(
            self, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the population additive genetic variance

        Parameters
        ----------
        Z : numpy.ndarray
            A matrix of genotypes.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population additive genetic 
            variances.

            Where:

            - ``t`` is the number of traits.
        """
        # estimate breeding values (n,t)
        gebv = self.gebv_numpy(Z, **kwargs)

        # calculate variance
        # (n,t).var(0) -> (t,)
        out = gebv.var(0)

        return out

    def var_A(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the population additive genetic variance

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population additive genetic 
            variances.

            Where:

            - ``t`` is the number of traits.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
        elif isinstance(gtobj, numpy.ndarray):
            Z = gtobj
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        out = self.var_A_numpy(Z, **kwargs)

        return out

    def var_a_numpy(
            self, 
            p: numpy.ndarray, 
            ploidy: Integral = 2, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : Integral
            Ploidy of the species.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population additive genic 
            variances.

            Where:

            - ``t`` is the number of traits.
        """
        # change shape to (p,1)
        p = p[:,None]

        # calculate additive genic variance
        # (p,t)**2 * (p,1) * (p,1) -> (p,t)
        # (p,t).sum[0] -> (t,)
        # scalar * (t,) -> (t,)
        out = (ploidy**2.0) * ((self.u_a**2) * p * (1.0 - p)).sum(0)

        return out

    def var_a(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            ploidy: Integral = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        ploidy : Integral, None
            Ploidy of the species.

            If ploidy is None:

            - If gtobj is a GenotypeMatrix, then get ploidy from GenotypeMatrix.
            - If gtobj is a numpy.ndarray, then assumed to be 2 (diploid).
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population additive genic 
            variances.

            Where:

            - ``t`` is the number of traits.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif isinstance(gtobj, numpy.ndarray):
            if ploidy is None:
                ploidy = 2
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate genic variance
        out = self.var_a_numpy(p, ploidy, **kwargs)

        return out

    def bulmer_numpy(
            self, 
            Z: numpy.ndarray, 
            p: numpy.ndarray, 
            ploidy: Integral = 2, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population Bulmer effect statistics.
            In the event that additive genic variance is zero, NaN's are produced.

            Where:

            - ``t`` is the number of traits.
        """
        sigma_A = self.var_A_numpy(Z)           # calculate additive genetic variance
        sigma_a = self.var_a_numpy(p, ploidy)   # calculate additive genetic variance
        mask = (sigma_a == 0.0)                 # determine where division by zero occurs
        denom = sigma_a.copy()                  # copy array
        denom[mask] = 1.0                       # substitute non-zero value
        out = sigma_A / denom                   # calculate Bulmer effect
        out[mask] = numpy.nan                   # add NaN's (avoids div by zero warning)
        return out

    def bulmer(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            ploidy: Optional[Integral] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
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

            - If gtobj is a GenotypeMatrix, then get ploidy from GenotypeMatrix.
            - If gtobj is a numpy.ndarray, then assumed to be 2 (diploid).
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population Bulmer effect statistics.
            In the event that additive genic variance is zero, NaN's are produced.

            Where:

            - ``t`` is the number of traits.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            Z = gtobj.mat_asformat("{0,1,2}")
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif isinstance(gtobj, numpy.ndarray):
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
    def usl_numpy(
            self, 
            p: numpy.ndarray, 
            ploidy: Integral, 
            unscale: bool = False, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the upper selection limit for a population.

        Parameters
        ----------
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : int
            Ploidy of the species.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population upper selection 
            limit statistics.

            Where:

            - ``t`` is the number of traits.
        """
        # reshape allele frequencies
        # (p,) -> (p,1)
        p = p[:,None]

        # get maximum attainable genotype
        # (p,t) ? (p,1) : (p,1) -> (p,t)
        uslgeno = numpy.where(
            self.u_a > 0.0,       # if the allele effect is positive
            p > 0.0,            # +allele: 1 if we have at least one +allele
            p >= 1.0            # -allele: 1 if we have fixation for -allele
        )

        # calculate usl value
        # scalar * (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = (float(ploidy) * self.u_a * uslgeno).sum(0)

        # make unscale adjustments if desired
        if unscale:
            # construct contrast (X*) to take cell means assuming that the first
            # column in X is a corner (first fixed effect is the intercept)

            # determine the number of fixed effects: 'q'
            nfixed = self.beta.shape[0]

            # allocate contrast matrix (X*)
            # (1,q)
            Xstar = numpy.empty(            # empty matrix
                (1, nfixed),                # (1,q)
                dtype = self.beta.dtype     # same dtype as beta
            )

            # fill first column with 1 since we assume it is the intercept
            Xstar[0,0] = 1

            # fill remaining columns with 1/q to calculate mean across all
            # fixed effects
            Xstar[0,1:] = 1 / nfixed

            # calculate intercepts (location)
            # (1,q) @ (q,t) -> (1,t)
            location = Xstar @ self.beta

            # add location to usl
            # (1,t) --ravel--> (t,)
            # (t,) + (t,) -> (t,)
            out += location.ravel()

        return out

    def usl(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            ploidy: Optional[Integral] = None, 
            unscale: bool = False, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the upper selection limit for a population.

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population upper selection 
            limit statistics.

            Where:

            - ``t`` is the number of traits.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif isinstance(gtobj, numpy.ndarray):
            if ploidy is None:
                ploidy = 2
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate genic variance
        out = self.usl_numpy(p, ploidy, unscale, **kwargs)

        return out

    def lsl_numpy(
            self, 
            p: numpy.ndarray, 
            ploidy: Integral, 
            unscale: bool = False, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the lower selection limit for a population.

        Parameters
        ----------
        p : numpy.ndarray
            A vector of genotype allele frequencies of shape (p,).
        ploidy : int
            Ploidy of the species.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population lower selection 
            limit statistics.

            Where:

            - ``t`` is the number of traits.
        """
        # reshape allele frequencies
        # (p,) -> (p,1)
        p = p[:,None]

        # get minimum attainable genotype
        # (p,t) ? (p,1) : (p,1) -> (p,t)
        lslgeno = numpy.where(
            self.u_a > 0.0,       # if the allele effect is positive
            p >= 1.0,           # +allele: 1 if we have fixation for +allele
            p > 0.0             # -allele: 1 if we have at least one -allele
        )

        # calculate lsl value
        # scalar * (p,t) * (p,t) -> (p,t)
        # (p,t).sum[0] -> (t,)
        out = (float(ploidy) * self.u_a * lslgeno).sum(0)

        # make unscale adjustments if desired
        if unscale:
            # construct contrast (X*) to take cell means assuming that the first
            # column in X is a corner (first fixed effect is the intercept)

            # determine the number of fixed effects: 'q'
            nfixed = self.beta.shape[0]

            # allocate contrast matrix (X*)
            # (1,q)
            Xstar = numpy.empty(            # empty matrix
                (1, nfixed),                # (1,q)
                dtype = self.beta.dtype     # same dtype as beta
            )

            # fill first column with 1 since we assume it is the intercept
            Xstar[0,0] = 1

            # fill remaining columns with 1/q to calculate mean across all
            # fixed effects
            Xstar[0,1:] = 1 / nfixed

            # calculate intercepts (location)
            # (1,q) @ (q,t) -> (1,t)
            location = Xstar @ self.beta

            # add location to usl
            # (1,t) --ravel--> (t,)
            # (t,) + (t,) -> (t,)
            out += location.ravel()

        return out

    def lsl(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            ploidy: Optional[Integral] = None, 
            unscale: bool = False, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the lower selection limit for a population.

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
            An array of shape ``(t,)`` containing population lower selection 
            limit statistics.

            Where:

            - ``t`` is the number of traits.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            p = gtobj.afreq()
            ploidy = gtobj.ploidy
        elif isinstance(gtobj, numpy.ndarray):
            if ploidy is None:
                ploidy = 2
            p = (1.0 / (ploidy * gtobj.shape[0])) * gtobj.sum(0)    # get allele frequencies
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        # calculate genic variance
        out = self.lsl_numpy(p, ploidy, unscale, **kwargs)

        return out

    ############ methods for allele attributes #############
    def facount(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Favorable allele count across all taxa.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to count favorable alleles.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native type.
        kwargs : dict
            Additional keyword arguments.
            
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,t)`` containing allele counts of the 
            favorable allele.
        """
        # process dtype
        if dtype is None:
            dtype = int
        dtype = numpy.dtype(dtype)

        # construct mask for beneficial alleles
        # (p,t)
        mask = (self.u_a > 0.0)

        # get allele counts for the genotype matrix
        # (p,) -> (p,1)
        acount = gmat.acount(dtype = dtype)[:,None]

        # get maximum number of favorable alleles
        # scalar
        maxfav = dtype.type(gmat.ploidy * gmat.ntaxa)

        # calculate favorable allele counts
        # (p,t)
        out = numpy.where(mask, acount, maxfav - acount)

        return out

    def fafreq(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Favorable allele frequency across all taxa.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to determine favorable allele frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native float 
            type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,t)`` containing allele frequencies of 
            the favorable allele.
        """
        # process dtype
        if dtype is None:
            dtype = float
        dtype = numpy.dtype(dtype)

        # get favorable allele frequencies
        out = (1.0 / (gmat.ploidy * gmat.ntaxa)) * self.facount(gmat)

        # convert datatypes if needed
        if out.dtype != dtype:
            out = out.astype(dtype)
        
        return out

    def faavail(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Determine which favorable alleles are available in an input set of taxa.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to determine favorable allele frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native bool type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,t)`` containing whether a favorable allele is available.
        """
        # process dtype
        if dtype is None:
            dtype = bool
        dtype = numpy.dtype(dtype)

        # get favorable allele counts
        facount = self.facount(gmat)

        # get boolean mask of favorable alleles that are available
        out = (facount > 0)

        # convert datatype if needed
        if out.dtype != dtype:
            out = out.astype(dtype)
        
        return out

    def fafixed(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Determine whether a favorable allele is fixed across all taxa.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to determine favorable allele frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,t)`` containing whether a favorable allele is fixed.
        """
        # process dtype
        if dtype is None:
            dtype = bool
        dtype = numpy.dtype(dtype)

        # get favorable allele counts
        facount = self.facount(gmat)

        # get maximum number of favorable alleles
        maxfav = gmat.ploidy * gmat.ntaxa

        # get boolean mask of favorable alleles that are fixed
        out = (facount == maxfav)

        # convert datatype if needed
        if out.dtype != dtype:
            out = out.astype(dtype)
        
        return out

    def dacount(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Deleterious allele count across all taxa.

        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to count deleterious alleles.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native type.
        kwargs : dict
            Additional keyword arguments.
            
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele counts of the deleterious allele.
        """
        # convert data type
        if dtype is None:
            dtype = int
        dtype = numpy.dtype(dtype)

        # construct mask for deleterious alleles
        # (p,t)
        mask = (self.u_a < 0.0)

        # get allele counts for the genotype matrix
        # (p,) -> (p,1)
        acount = gmat.acount(dtype = dtype)[:,None]

        # get maximum number of favorable alleles
        # scalar
        maxfav = dtype.type(gmat.ploidy * gmat.ntaxa)

        # calculate favorable allele counts
        # (p,t)
        out = numpy.where(mask, acount, maxfav - acount)

        return out

    def dafreq(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Deleterious allele frequency across all taxa.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to determine deleterious allele frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies of the deleterious allele.
        """
        # process dtype
        if dtype is None:
            dtype = float
        dtype = numpy.dtype(dtype)

        # get deleterious allele frequencies
        out = (1.0 / (gmat.ploidy * gmat.ntaxa)) * self.dacount(gmat)

        # convert datatypes if needed
        if out.dtype != dtype:
            out = out.astype(dtype)

        return out
    
    def daavail(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Determine whether a deleterious allele is available in the present taxa.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to determine deleterious allele frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native boolean type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,t)`` containing whether a deleterious allele is available.
        """
        # process dtype
        if dtype is None:
            dtype = bool
        dtype = numpy.dtype(dtype)

        # get deleterious allele counts
        facount = self.dacount(gmat)

        # get boolean mask of deleterious alleles that are available
        out = (facount > 0)

        # convert datatype if needed
        if out.dtype != dtype:
            out = out.astype(dtype)
        
        return out

    def dafixed(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Determine whether a deleterious allele is fixed across all taxa.
        
        Parameters
        ----------
        gmat : GenotypeMatrix
            Genotype matrix for which to determine deleterious allele frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,)`` containing whether a deleterious allele is fixed.
        """
        # process dtype
        if dtype is None:
            dtype = bool
        dtype = numpy.dtype(dtype)

        # get deleterious allele counts
        dacount = self.dacount(gmat)

        # get maximum number of deleterious alleles
        maxdel = gmat.ploidy * gmat.ntaxa

        # get boolean mask of deleterious alleles that are fixed
        out = (dacount == maxdel)

        # convert datatype if needed
        if out.dtype != dtype:
            out = out.astype(dtype)
        
        return out

    ################## Model I/O methods ###################
    def to_pandas_dict(
            self,
            trait_cols: Optional[Union[str,Sequence]] = "trait",
            **kwargs: dict
        ) -> Dict[str,pandas.DataFrame]:
        """
        Export a DenseAdditiveLinearGenomicModel to a ``dict`` of ``pandas.DataFrame``.

        Parameters
        ----------
        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to write regression coefficients.
            If ``Sequence``, column names are given by the strings in the 
            ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"trait"``. Use trait names given in 
            the ``trait`` property.
            If ``None``, use numeric trait column names.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a 
            ``dict`` of ``pandas.DataFrame``.

        Returns
        -------
        out : dict
            An output dataframe.
        """
        # type checks
        if trait_cols is not None:
            if isinstance(trait_cols, str):
                check_str_value(trait_cols, "trait_cols", "trait")
            elif isinstance(trait_cols, Sequence):
                check_len(trait_cols, "trait_cols", self.ntrait)
            else:
                check_is_str_or_Sequence(trait_cols, "trait_cols")

        # process trait_cols
        if trait_cols is None:
            trait_cols = numpy.arange(self.ntrait)
        elif isinstance(trait_cols, str):
            trait_cols = numpy.arange(self.ntrait) if self.trait is None else self.trait

        # construct output dictionary
        out = {
            "beta": pandas.DataFrame(self.beta, columns = trait_cols),
            "u_misc": pandas.DataFrame(self.u_misc, columns = trait_cols),
            "u_a": pandas.DataFrame(self.u_a, columns = trait_cols),
        }

        return out

    def to_csv_dict(
            self,
            filenames: Dict[str,str],
            trait_cols: Optional[Union[str,Sequence]] = "trait",
            sep: str = ',', 
            header: bool = True, 
            index: bool = False, 
            **kwargs: dict
        ) -> None:
        """
        Export a DenseAdditiveLinearGenomicModel to a set of CSV files specified 
        by values in a ``dict``.

        Parameters
        ----------
        filenames : dict of str
            CSV file names to which to write. Must have the keys: ``"beta"``, 
            ``"u_misc"``, and ``"u_a"`` (case sensitive).

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to write regression coefficients.
            If ``Sequence``, column names are given by the strings in the 
            ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"trait"``. Use trait names given in 
            the ``trait`` property.
            If ``None``, use numeric trait column names.

        sep : str, default = ","
            Separator to use in the exported CSV files.
        
        header : bool, default = True
            Whether to save header names.
        
        index : bool, default = False
            Whether to save a row index in the exported CSV files.

        kwargs : dict
            Additional keyword arguments to use for dictating export to a CSV.
        """
        # type checks
        check_is_dict(filenames, "filenames")
        check_dict_has_keys(filenames, "filenames", "beta", "u_misc", "u_a")

        # export to dict of pandas.DataFrame
        df_dict = self.to_pandas_dict(
            trait_cols = trait_cols,
        )

        # export to CSV using pandas
        if filenames["beta"] is not None:
            df_dict["beta"].to_csv(
                path_or_buf = filenames["beta"],
                sep = sep,
                header = header,
                index = index,
                **kwargs
            )
        if filenames["u_misc"] is not None:
            df_dict["u_misc"].to_csv(
                path_or_buf = filenames["u_misc"],
                sep = sep,
                header = header,
                index = index,
                **kwargs
            )
        if filenames["u_a"] is not None:
            df_dict["u_a"].to_csv(
                path_or_buf = filenames["u_a"],
                sep = sep,
                header = header,
                index = index,
                **kwargs
            )

    def to_hdf5(
            self, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> None:
        """
        Write GenomicModel to an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name to which to write.
        groupname : str or None
            HDF5 group name under which GenomicModel data is stored.
            If ``None``, GenomicModel is written to the base HDF5 group.
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
            "u_misc": self.u_misc,
            "u_a": self.u_a,
            "trait": self.trait,
            "model_name": self.model_name,
            "hyperparams": self.hyperparams
        }
        save_dict_to_hdf5(h5file, groupname, data_dict)         # write data
        ######################################################### write conclusion
        h5file.close()                                          # close the file

    ############################## Class Methods ###############################

    ################## Model I/O methods ###################
    @classmethod
    def from_pandas_dict(
            cls,
            dic: Dict[str,pandas.DataFrame],
            trait_cols: Optional[Union[str,Sequence]] = "infer",
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Read an object from a ``dict`` of ``pandas.DataFrame``.

        Parameters
        ----------
        dic : dict
            Python dictionary containing ``pandas.DataFrame`` from which to read.
            Must have the following fields::

            - ``"beta"`` - ``pandas.DataFrame`` containing fixed effects.
            - ``"u_misc"`` - ``None`` or ``pandas.DataFrame`` containing 
              miscellaneous random effects.
            - ``"u_a"`` - ``None`` or ``pandas.DataFrame`` containing additive 
              genetic marker random effects.

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to read regression coefficients.
            If ``Sequence``, column names are given by the strings or integers 
            in the ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"infer"``. Use columns in the 
            ``"beta"`` input dataframe to load trait breeding values.
            If ``None``, do not load any trait regression coefficients.

        model_name : str, None
            Name of the model.

        hyperparams : dict, None
            Model parameters.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a 
            ``dict`` of ``pandas.DataFrame``.

        Returns
        -------
        out : DenseAdditiveLinearGenomicModel
            A DenseAdditiveLinearGenomicModel read from a ``dict`` of ``pandas.DataFrame``.
        """
        # type checks
        check_is_dict(dic, "dic")
        check_dict_has_keys(dic, "dic", "beta", "u_misc", "u_a")
        check_is_pandas_DataFrame(dic["beta"], 'dic["beta"]')
        if dic["u_misc"] is not None:
            check_is_pandas_DataFrame(dic["u_misc"], 'dic["u_misc"]')
        if dic["u_a"] is not None:
            check_is_pandas_DataFrame(dic["u_a"], 'dic["u_a"]')
        if trait_cols is not None:
            if isinstance(trait_cols, str):
                check_str_value(trait_cols, "trait_cols", "infer")
            elif isinstance(trait_cols, Sequence):
                pass
            else:
                check_is_str_or_Sequence(trait_cols, "trait_cols")

        # get trait values from dic["beta"]
        trait = None
        if trait_cols is None:
            pass
        elif isinstance(trait_cols, str):
            if trait_cols == "infer":
                trait = dic["beta"].columns.to_numpy(dtype = object)
        elif isinstance(trait_cols, Sequence):
            trait = numpy.array(trait_cols, dtype = object)
        
        # get beta values 
        df = dic["beta"]
        ix = [] if trait is None else [df.columns.get_loc(e) if isinstance(e,str) else e for e in trait]
        beta = df.iloc[:,ix].to_numpy(dtype = float)
        
        # get u_misc values
        u_misc = None
        if dic["u_misc"] is not None:
            df = dic["u_misc"]
            ix = [] if trait is None else [df.columns.get_loc(e) if isinstance(e,str) else e for e in trait]
            u_misc = df.iloc[:,ix].to_numpy(dtype = float)
        
        # get u_a values
        u_a = None
        if dic["u_a"] is not None:
            df = dic["u_a"]
            ix = [] if trait is None else [df.columns.get_loc(e) if isinstance(e,str) else e for e in trait]
            u_a = df.iloc[:,ix].to_numpy(dtype = float)

        # create output object
        out = cls(
            beta = beta,
            u_misc = u_misc,
            u_a = u_a,
            trait = trait,
            model_name = model_name,
            hyperparams = hyperparams,
        )

        return out

    @classmethod
    def from_csv_dict(
            cls, 
            filenames: Dict[str,str],
            sep: str = ',',
            header: int = 0,
            trait_cols: Optional[Union[str,Sequence]] = "infer",
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Read a DenseAdditiveLinearGenomicModel from a set of CSV files specified by values in a ``dict``.

        Parameters
        ----------
        filename : str
            Dictionary of CSV file names from which to read.
            Must have the following fields::

            - ``"beta"`` - ``str`` containing fixed effects.
            - ``"u_misc"`` - ``None`` or ``str`` of CSV file path containing 
              miscellaneous random effects.
            - ``"u_a"`` - ``None`` or ``str`` of CSV file path containing additive 
              genetic marker random effects.

        sep : str, default = ','
            CSV delimiter to use.
        
        header : int, list of int, default=0
            Row number(s) to use as the column names, and the start of the data.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        trait_cols : Sequence, str, None, default = "trait"
            Names of the trait columns to which to read regression coefficients.
            If ``Sequence``, column names are given by the strings or integers 
            in the ``trait_cols`` Sequence.
            If ``str``, must be equal to ``"infer"``. Use columns in the 
            ``"beta"`` input dataframe to load trait breeding values.
            If ``None``, do not load any trait regression coefficients.

        model_name : str, None
            Name of the model.

        hyperparams : dict, None
            Model parameters.

        kwargs : dict
            Additional keyword arguments to use for dictating importing from a CSV.

        Returns
        -------
        out : DenseAdditiveLinearGenomicModel
            A DenseAdditiveLinearGenomicModel read from a set of CSV files.
        """
        # type checks
        check_is_dict(filenames, "filenames")
        check_dict_has_keys(filenames, "filenames", "beta", "u_misc", "u_a")
        check_is_str(filenames["beta"], 'filenames["beta"]')
        if filenames["u_misc"] is not None:
            check_is_str(filenames["u_misc"], 'filenames["u_misc"]')
        if filenames["u_a"] is not None:
            check_is_str(filenames["u_a"], 'filenames["u_a"]')

        # read files
        beta_df = pandas.read_csv(filenames["beta"], sep = sep, header = header, **kwargs)
        u_misc_df = None if filenames["u_misc"] is None else pandas.read_csv(filenames["u_misc"], sep = sep, header = header, **kwargs)
        u_a_df = None if filenames["u_a"] is None else pandas.read_csv(filenames["u_a"], sep = sep, header = header, **kwargs)

        # create pandas dictionary
        dic = {
            "beta": beta_df,
            "u_misc": u_misc_df,
            "u_a": u_a_df,
        }

        # create output
        out = cls.from_pandas_dict(
            dic = dic, 
            trait_cols = trait_cols, 
            model_name = model_name, 
            hyperparams = hyperparams, 
        )

        return out

    @classmethod
    def from_hdf5(
            cls, 
            filename: str, 
            groupname: Optional[str] = None
        ) -> 'DenseAdditiveLinearGenomicModel':
        """
        Read GenotypeMatrix from an HDF5 file.

        Parameters
        ----------
        filename : str
            HDF5 file name which to read.
        groupname : str or None
            HDF5 group name under which GenotypeMatrix data is stored.
            If ``None``, GenotypeMatrix is read from base HDF5 group.

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
        required_fields = ["beta", "u_misc", "u_a"]             # all required arguments
        for field in required_fields:                           # for each required field
            fieldname = groupname + field                       # concatenate base groupname and field
            check_group_in_hdf5(fieldname, h5file, filename)    # check that group exists
        ######################################################### read data
        data_dict = {                                           # output dictionary
            "beta": None,
            "u_misc" : None,
            "u_a" : None,
            "trait": None,
            "model_name": None,
            "hyperparams": None
        }
        data_dict["beta"] = h5file[groupname + "beta"][()]      # read beta array
        data_dict["u_misc"] = h5file[groupname + "u_misc"][()]  # read u array
        data_dict["u_a"] = h5file[groupname + "u_a"][()]        # read u array
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
        fieldname = groupname + "hyperparams"                        # construct "groupname/hyperparams"
        if fieldname in h5file:                                 # if "groupname/hyperparams" in hdf5
            data_dict["hyperparams"] = {}                            # create empty dictionary
            view = h5file[fieldname]                            # get view of dataset
            for key in view.keys():                             # for each field
                data_dict["hyperparams"][key] = view[key][()]        # extract data
        ######################################################### read conclusion
        h5file.close()                                          # close file
        ######################################################### create object
        dalgmod = cls(**data_dict)                              # create object from read data
        return dalgmod



################################## Utilities ###################################
def check_is_DenseAdditiveLinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type DenseAdditiveLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseAdditiveLinearGenomicModel):
        raise TypeError("variable '{0}' must be a DenseAdditiveLinearGenomicModel".format(vname))
