"""
Module implementign classes and error checking routines for genomic
prediciton models that incorporate additive and dominance effects.
"""

import copy
from pathlib import Path
import h5py
from numbers import Integral
from typing import Dict, Optional, Sequence, Union
import numpy
import pandas
from pybrops.core.error.error_io_python import check_file_exists
from pybrops.core.error.error_type_numpy import check_is_ndarray, check_ndarray_dtype_is_float64, check_ndarray_dtype_is_object
from pybrops.core.error.error_type_pandas import check_is_pandas_DataFrame
from pybrops.core.error.error_type_python import check_is_dict, check_is_str, check_is_str_or_Sequence
from pybrops.core.error.error_value_h5py import check_h5py_File_has_group, check_h5py_File_is_readable, check_h5py_File_is_writable
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_eq, check_ndarray_ndim
from pybrops.core.error.error_value_python import check_dict_has_keys, check_len, check_str_value
from pybrops.core.util.h5py import h5py_File_read_dict, h5py_File_read_ndarray, h5py_File_read_ndarray_utf8, h5py_File_read_utf8, h5py_File_write_dict

from pybrops.model.gmod.AdditiveDominanceLinearGenomicModel import AdditiveDominanceLinearGenomicModel
from pybrops.model.gmod.DenseAdditiveLinearGenomicModel import DenseAdditiveLinearGenomicModel
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.bvmat.DenseGenomicEstimatedBreedingValueMatrix import DenseGenomicEstimatedBreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix


class DenseAdditiveDominanceLinearGenomicModel(
        DenseAdditiveLinearGenomicModel,
        AdditiveDominanceLinearGenomicModel,
    ):
    """
    The DenseAdditiveDominanceLinearGenomicModel class represents a Multivariate Multiple
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
        \\mathbf{Z} = \\begin{bmatrix} \\mathbf{Z_{misc}} & \\mathbf{Z_{a}} & \\mathbf{Z_{d}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{Z_{misc}}` is a matrix of miscellaneous random effect predictors of shape ``(n,p_misc)``
    - :math:`\\mathbf{Z_{a}}` is a matrix of additive genomic marker predictors of shape ``(n,p_a)``
    - :math:`\\mathbf{Z_{d}}` is a matrix fo dominance genomic marker predictors of shape ``(n,p_d)``

    .. math::
        \\mathbf{U} = \\begin{bmatrix} \\mathbf{U_{misc}} \\\\ \\mathbf{U_{a}} \\\\
        \\mathbf{U_{d}} \\end{bmatrix}

    Where:

    - :math:`\\mathbf{U_{misc}}` is a matrix of miscellaneous random effects of shape ``(p_misc,t)``
    - :math:`\\mathbf{U_{a}}` is a matrix of additive genomic marker effects of shape ``(p_a,t)``
    - :math:`\\mathbf{U_{d}}` is a matrix of dominance genomic marker effects of shape ``(p_d,t)``

    Shape definitions:

    - ``n`` is the number of individuals
    - ``q`` is the number of fixed effect predictors (e.g. environments)
    - ``p`` is the number of random effect predictors.
    - ``p_misc`` is the number of miscellaneous random effect predictors.
    - ``p_a`` is the number of additive genomic marker predictors.
    - ``p_d`` is the number of dominance genomic marker predictors.
    - The sum of ``p_misc`` and ``p_a`` and ``p_d`` equals ``p``.
    - ``t`` is the number of traits
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            beta: numpy.ndarray, 
            u_misc: Union[numpy.ndarray,None], 
            u_a: Union[numpy.ndarray,None], 
            u_d: Union[numpy.ndarray,None], 
            trait: Optional[numpy.ndarray] = None, 
            model_name: Optional[str] = None, 
            hyperparams: Optional[dict] = None, 
            **kwargs: dict
        ) -> None:
        """
        Constructor for DenseAdditiveDominanceLinearGenomicModel class.

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
        
        u_d : numpy.ndarray, None
            A ``float64`` random effect regression coefficient matrix of shape ``(p_d,t)`` containing dominance marker effects.

            Where:

            - ``p_d`` is the number of dominance marker effect predictors. Must be equal to ``p_a``.
            - ``t`` is the number of traits.

            If ``None``, then set a zero array of shape ``(p_a,t)``
        
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
        # call super() constructor
        super(DenseAdditiveDominanceLinearGenomicModel, self).__init__(
            beta = beta,
            u_misc = u_misc,
            u_a = u_a,
            trait = trait,
            model_name = model_name,
            hyperparams = hyperparams,
            **kwargs
        )
        # set dominance factors
        self.u_d = u_d

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
        return "<{0} of shape [B = {1}, U_misc = {2}, U_a = {3}, U_d = {4}] at {5}>".format(
            type(self).__name__,
            self.beta.shape,
            self.u_misc.shape,
            self.u_a.shape,
            self.u_d.shape,
            hex(id(self))
        )

    def __copy__(
            self
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : DenseAdditiveDominanceLinearGenomicModel
            A shallow copy of the model.
        """
        out = self.__class__(
            beta        = copy.copy(self.beta),
            u_misc      = copy.copy(self.u_misc),
            u_a         = copy.copy(self.u_a),
            u_d         = copy.copy(self.u_d),
            trait       = copy.copy(self.trait),
            model_name  = copy.copy(self.model_name),
            hyperparams = copy.copy(self.hyperparams),
        )

        return out

    def __deepcopy__(
            self, 
            memo: Optional[dict] = None
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Make a deep copy of the GenomicModel.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseAdditiveDominanceLinearGenomicModel
            A deep copy of the model.
        """
        out = self.__class__(
            beta        = copy.deepcopy(self.beta, memo),
            u_misc      = copy.deepcopy(self.u_misc, memo),
            u_a         = copy.deepcopy(self.u_a, memo),
            u_d         = copy.deepcopy(self.u_d, memo),
            trait       = copy.deepcopy(self.trait, memo),
            model_name  = copy.deepcopy(self.model_name, memo),
            hyperparams = copy.deepcopy(self.hyperparams, memo),
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

    # inherit nexplan_beta
    # inherit nparam_beta
    # inherit beta

    @property
    def nexplan_u(self) -> Integral:
        """Number of random effect explanatory variables required by the model."""
        return self.nexplan_u_misc + self.nexplan_u_a + self.nexplan_u_d

    @property
    def nparam_u(self) -> Integral:
        """Number of random effect parameters."""
        return self.nparam_u_misc + self.nparam_u_a + self.nparam_u_d

    @property
    def u(self) -> numpy.ndarray:
        """Random effect regression coefficients."""
        # get random effects and concatenate along compatible axes
        out = numpy.concatenate([self.u_misc, self.u_a, self.u_d], axis = 0)
        return out
    @u.setter
    def u(self, value: numpy.ndarray) -> None:
        """Set random effect regression coefficients"""
        raise AttributeError("variable 'u' is read-only; use 'u_misc' and 'u_a' to modify 'u'.")

    ####### Additive Linear Genomic Model Parameters #######

    # inherit nexplan_u_misc
    # inherit nparam_u_misc
    # inherit u_misc

    # inherit nexplan_u_a
    # inherit nparam_u_a
    # inherit u_a

    @property
    def nexplan_u_d(self) -> Integral:
        """Number of dominance genomic marker explanatory variables required by the model."""
        return self.u_d.shape[0]

    @property
    def nparam_u_d(self) -> Integral:
        """Number of additive genomic marker parameters."""
        return self.u_d.size

    @property
    def u_d(self) -> numpy.ndarray:
        """Additive genomic marker effects."""
        return self._u_d
    @u_d.setter
    def u_d(self, value: numpy.ndarray) -> None:
        """Set additive genomic marker effect regression coefficients"""
        if value is None:                                   # if value is None
            t = self.beta.shape[1]                          # get number of traits
            p_d = self.u_a.shape[0]                         # get number of markers
            value = numpy.zeros((p_d,t), dtype = float)     # make zero array of shape (p_d,t)
        check_is_ndarray(value, "u_d")
        check_ndarray_ndim(value, "u_d", 2)
        check_ndarray_dtype_is_float64(value, "u_d")
        self._u_d = value

    ################## Genomic Model Data ##################
    
    # inherit model_name
    # inherit hyperparams
    # inherit trait
    # inherit ntrait

    ############################## Object Methods ##############################

    #################### Model copying #####################
    def copy(
            self
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Make a shallow copy of the DenseAdditiveDominanceLinearGenomicModel.

        Returns
        -------
        out : DenseAdditiveDominanceLinearGenomicModel
            A shallow copy of the original DenseAdditiveDominanceLinearGenomicModel
        """
        return self.__copy__()

    def deepcopy(
            self,
            memo: Optional[dict] = None
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Make a deep copy of the DenseAdditiveDominanceLinearGenomicModel.

        Parameters
        ----------
        memo : dict
            Dictionary of memo metadata.

        Returns
        -------
        out : DenseAdditiveDominanceLinearGenomicModel
            A deep copy of the original DenseAdditiveDominanceLinearGenomicModel
        """
        return self.__deepcopy__(memo)

    ####### methods for model fitting and prediction #######
    @classmethod
    def fit_numpy(
            cls, 
            Y: numpy.ndarray, 
            X: numpy.ndarray, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
        """
        Fit a dense, additive + dominance linear genomic model.

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
        raise AttributeError("DenseAdditiveDominanceLinearGenomicModel is read-only")

    @classmethod
    def fit(
            cls, 
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
        raise AttributeError("DenseAdditiveDominanceLinearGenomicModel is read-only")

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
            A matrix of genotype values coded as {0,1,2} for additive predictors and {0,1} for dominance predictors.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        Y_hat : numpy.ndarray
            A matrix of estimated breeding values.
        """
        # input checks
        check_is_ndarray(X, "X")
        check_ndarray_ndim(X, "X", 2)
        check_ndarray_axis_len_eq(X, "X", 1, self.nexplan_beta)
        check_is_ndarray(Z, "Z")
        check_ndarray_ndim(Z, "Z", 2)
        check_ndarray_axis_len_eq(Z, "Z", 1, self.nexplan_u)

        # Y = XB + ZU
        # (n,t) = (n,q) @ (q,t) + (n,p) @ (p,t)
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
            values. If ``numpy.ndarray``, must be coded as ``{0,1,2}``.
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
            A = gtobj.mat_asformat("{0,1,2}")                   # get allele counts
            D = numpy.logical_and(A != 0, A != gtobj.ploidy)    # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
            taxa = gtobj.taxa
            taxa_grp = gtobj.taxa_grp
        elif isinstance(gtobj, numpy.ndarray):
            A = gtobj                                           # get allele counts
            D = gtobj == 1                                      # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
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
        Y_hat = self.predict_numpy(X, Z, **kwargs)

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
        elif isinstance(ptobj, numpy.ndarray):
            Y = ptobj
        else:
            raise TypeError("must be BreedingValueMatrix, numpy.ndarray")

        # process cvobj
        if isinstance(cvobj, numpy.ndarray):
            X = cvobj
        else:
            raise TypeError("must be numpy.ndarray")

        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            A = gtobj.mat_asformat("{0,1,2}")                   # get allele counts
            D = numpy.logical_and(A != 0, A != gtobj.ploidy)    # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
        elif isinstance(gtobj, numpy.ndarray):
            A = gtobj                                           # get allele counts
            D = gtobj == 1                                      # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
        else:
            raise TypeError("must be GenotypeMatrix, numpy.ndarray")

        # calculate coefficient of determination
        Rsq = self.score_numpy(Y, X, Z, **kwargs)

        return Rsq

    ######## methods for estimated breeding values #########
    
    # inherit gebv_numpy
    # inherit gebv

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
        # type checks
        check_is_ndarray(Z, "Z")
        check_ndarray_ndim(Z, "Z", 2)
        check_ndarray_axis_len_eq(Z, "Z", 1, self.nexplan_u_a + self.nexplan_u_d)

        # concatenate U
        U = numpy.concatenate([self.u_a, self.u_d], axis = 0)

        # Y = ZU
        gegv_hat = (Z @ U)

        return gegv_hat
    
    def gegv(
            self,
            gtobj: Union[GenotypeMatrix,numpy.ndarray],
            **kwargs: dict
        ) -> BreedingValueMatrix:
        """
        Calculate genomic estimated genotypic values.

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            A matrix of genotypic markers.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A matrix of genomic estimated genotypic values.
        """
        # process gtobj
        if isinstance(gtobj, GenotypeMatrix):
            A = gtobj.mat_asformat("{0,1,2}")                   # get allele counts
            D = numpy.logical_and(A != 0, A != gtobj.ploidy)    # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
            taxa = gtobj.taxa
            taxa_grp = gtobj.taxa_grp
        elif isinstance(gtobj, numpy.ndarray):
            A = gtobj                                           # get allele counts
            D = gtobj == 1                                      # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
            taxa = None
            taxa_grp = None
        else:
            raise TypeError("accepted types are GenotypeMatrix, numpy.ndarray")

        # make predictions
        gegv_hat = self.gegv_numpy(Z, **kwargs)

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
        gegv_hat += location

        # create output breeding value matrix
        out = DenseGenomicEstimatedBreedingValueMatrix.from_numpy(
            mat = gegv_hat,
            taxa = taxa,
            taxa_grp = taxa_grp,
            trait = self.trait
        )

        return out

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
        gegv = self.gegv_numpy(Z, **kwargs)

        # calculate variance
        # (n,t).var(0) -> (t,)
        out = gegv.var(0)

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
            A = gtobj.mat_asformat("{0,1,2}")                   # get allele counts
            D = numpy.logical_and(A != 0, A != gtobj.ploidy)    # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
        elif isinstance(gtobj, numpy.ndarray):
            A = gtobj                                           # get allele counts
            D = gtobj == 1                                      # get heterozygous loci indicator variables
            Z = numpy.concatenate([A, D], axis = 1)             # create additive + dominance
        else:
            raise TypeError("must be GenotypeMatrix, ndarray")

        out = self.var_G_numpy(Z, **kwargs)

        return out

    # inherit var_A_numpy
    # inherit var_A
    # inherit var_a_numpy
    # inherit var_a
    # inherit bulmer_numpy
    # inherit bulmer

    ############# methods for selection limits #############

    # inherit usl_numpy
    # inherit usl
    # inherit lsl_numpy
    # inherit lsl

    ############ methods for allele attributes #############

    # inherit facount
    # inherit fafreq
    # inherit faavail
    # inherit fafixed
    # inherit dacount
    # inherit dafreq
    # inherit daavail
    # inherit dafixed

    ################## Model I/O methods ###################
    def to_pandas_dict(
            self,
            trait_cols: Optional[Union[str,Sequence]] = "trait",
            **kwargs: dict
        ) -> Dict[str,pandas.DataFrame]:
        """
        Export a DenseAdditiveDominanceLinearGenomicModel to a ``dict`` of ``pandas.DataFrame``.

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
            "beta":   pandas.DataFrame(self.beta,   columns = trait_cols),
            "u_misc": pandas.DataFrame(self.u_misc, columns = trait_cols),
            "u_a":    pandas.DataFrame(self.u_a,    columns = trait_cols),
            "u_d":    pandas.DataFrame(self.u_d,    columns = trait_cols),
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
        Export a DenseAdditiveDominanceLinearGenomicModel to a set of CSV files specified 
        by values in a ``dict``.

        Parameters
        ----------
        filenames : dict of str
            CSV file names to which to write. Must have the keys: ``"beta"``, 
            ``"u_misc"``, ``"u_a"``, and ``"u_d"`` (case sensitive).

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
        check_dict_has_keys(filenames, "filenames", "beta", "u_misc", "u_a", "u_d")

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
        if filenames["u_d"] is not None:
            df_dict["u_d"].to_csv(
                path_or_buf = filenames["u_d"],
                sep = sep,
                header = header,
                index = index,
                **kwargs
            )

    def to_hdf5(
            self, 
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None,
            overwrite: bool = True,
        ) -> None:
        """
        Write ``DenseAdditiveDominanceLinearGenomicModel`` to an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str``, an HDF5 file name to which to write. File is closed after writing.
            If ``h5py.File``, an opened HDF5 file to which to write. File is not closed after writing.

        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseAdditiveDominanceLinearGenomicModel`` data is stored.
            If ``None``, ``DenseAdditiveDominanceLinearGenomicModel`` is written to the base HDF5 group.

        overwrite : bool
            Whether to overwrite data fields if they are present in the HDF5 file.
        """
        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r+``) mode
        if isinstance(filename, (str,Path)):
            h5file = h5py.File(filename, "a")

        # elif we have an h5py.File, make sure mode is writable, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_writable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        #### write data to HDF5 file and (optionally) close ####

        # data dictionary
        data = {
            "beta"          : self.beta,
            "u_misc"        : self.u_misc,
            "u_a"           : self.u_a,
            "u_d"           : self.u_d,
            "trait"         : self.trait,
            "model_name"    : self.model_name,
            "hyperparams"   : self.hyperparams,
        }

        # save data
        h5py_File_write_dict(h5file, groupname, data, overwrite)

        # close the file, only if the provided filename was a string or Path and not a h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

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
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Read an object from a ``dict`` of ``pandas.DataFrame``.

        Parameters
        ----------
        dic : dict
            Python dictionary containing ``pandas.DataFrame`` from which to read.
            Must have the following fields::

            - ``"beta"`` is a ``pandas.DataFrame`` containing fixed effects.
            - ``"u_misc"`` is ``None`` or ``pandas.DataFrame`` containing 
              miscellaneous random effects.
            - ``"u_a"`` is ``None`` or a ``pandas.DataFrame`` containing additive 
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
        out : DenseAdditiveDominanceLinearGenomicModel
            A DenseAdditiveDominanceLinearGenomicModel read from a ``dict`` of ``pandas.DataFrame``.
        """
        # type checks
        check_is_dict(dic, "dic")
        check_dict_has_keys(dic, "dic", "beta", "u_misc", "u_a", "u_d")
        check_is_pandas_DataFrame(dic["beta"], 'dic["beta"]')
        if dic["u_misc"] is not None:
            check_is_pandas_DataFrame(dic["u_misc"], 'dic["u_misc"]')
        if dic["u_a"] is not None:
            check_is_pandas_DataFrame(dic["u_a"], 'dic["u_a"]')
        if dic["u_d"] is not None:
            check_is_pandas_DataFrame(dic["u_d"], 'dic["u_d"]')
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

        # get u_d values
        u_d = None
        if dic["u_d"] is not None:
            df = dic["u_d"]
            ix = [] if trait is None else [df.columns.get_loc(e) if isinstance(e,str) else e for e in trait]
            u_d = df.iloc[:,ix].to_numpy(dtype = float)

        # create output object
        out = cls(
            beta        = beta,
            u_misc      = u_misc,
            u_a         = u_a,
            u_d         = u_d,
            trait       = trait,
            model_name  = model_name,
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
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Read a ``DenseAdditiveDominanceLinearGenomicModel`` from a set of CSV files specified by values in a ``dict``.

        Parameters
        ----------
        filenames : str
            Dictionary of CSV file names from which to read.
            
            Must have the following fields::

            - ``"beta"`` is a ``str`` containing fixed effects.
            - ``"u_misc"`` is ``None`` or a ``str`` of CSV file path containing 
              miscellaneous random effects.
            - ``"u_a"`` is ``None`` or a ``str`` of CSV file path containing additive 
              genetic marker random effects.
            - ``"u_d"`` is ``None`` or a ``str`` of CSV file path containing dominance 
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
        out : DenseAdditiveDominanceLinearGenomicModel
            A ``DenseAdditiveDominanceLinearGenomicModel`` read from a set of CSV files.
        """
        # type checks
        check_is_dict(filenames, "filenames")
        check_dict_has_keys(filenames, "filenames", "beta", "u_misc", "u_a", "u_d")
        check_is_str(filenames["beta"], 'filenames["beta"]')
        if filenames["u_misc"] is not None:
            check_is_str(filenames["u_misc"], 'filenames["u_misc"]')
        if filenames["u_a"] is not None:
            check_is_str(filenames["u_a"], 'filenames["u_a"]')
        if filenames["u_d"] is not None:
            check_is_str(filenames["u_d"], 'filenames["u_d"]')

        # read files
        beta_df = pandas.read_csv(filenames["beta"], sep = sep, header = header, **kwargs)
        u_misc_df = None if filenames["u_misc"] is None else pandas.read_csv(filenames["u_misc"], sep = sep, header = header, **kwargs)
        u_a_df = None if filenames["u_a"] is None else pandas.read_csv(filenames["u_a"], sep = sep, header = header, **kwargs)
        u_d_df = None if filenames["u_d"] is None else pandas.read_csv(filenames["u_d"], sep = sep, header = header, **kwargs)

        # create pandas dictionary
        dic = {
            "beta":   beta_df,
            "u_misc": u_misc_df,
            "u_a":    u_a_df,
            "u_d":    u_d_df,
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
            filename: Union[str,Path,h5py.File], 
            groupname: Optional[str] = None
        ) -> 'DenseAdditiveDominanceLinearGenomicModel':
        """
        Read ``DenseAdditiveDominanceLinearGenomicModel`` from an HDF5 file.

        Parameters
        ----------
        filename : str, Path, h5py.File
            If ``str`` or ``Path``, an HDF5 file name from which to read. File is closed after reading.
            If ``h5py.File``, an opened HDF5 file from which to read. File is not closed after reading.
        groupname : str, None
            If ``str``, an HDF5 group name under which ``DenseAdditiveDominanceLinearGenomicModel`` data is stored.
            If ``None``, ``DenseAdditiveDominanceLinearGenomicModel`` is read from base HDF5 group.

        Returns
        -------
        out : DenseAdditiveDominanceLinearGenomicModel
            A genomic model read from file.
        """
        ########################################################
        ############ process ``groupname`` argument ############

        # if we have a string
        if isinstance(groupname, str):
            # if last character in string is not '/', add '/' to end of string
            if groupname[-1] != '/':
                groupname += '/'
        
        # else if ``groupname`` is None, set ``groupname`` to empty string
        elif groupname is None:
            groupname = ""
        
        # else raise error
        else:
            raise TypeError(
                "``groupname`` must be of type ``str`` or ``None`` but received type ``{0}``".format(
                    type(groupname).__name__
                )
            )

        ########################################################
        ############ process ``filename`` argument #############

        # HDF5 file object
        h5file = None

        # if we have a string or Path, open HDF5 file in append (``r``) mode
        if isinstance(filename, (str,Path)):
            check_file_exists(filename)
            h5file = h5py.File(filename, "r")

        # elif we have an h5py.File, make sure mode is in at least ``r`` mode, and copy pointer
        elif isinstance(filename, h5py.File):
            check_h5py_File_is_readable(filename)
            h5file = filename
        
        # else raise TypeError
        else:
            raise TypeError(
                "``filename`` must be of type ``str``, ``Path``, or ``h5py.File`` but received type ``{0}``".format(
                    type(filename).__name__
                )
            )

        ########################################################
        ######## check that we have all required fields ########

        # all required arguments
        required_fields = ["beta", "u_misc", "u_a", "u_d"]

        # for each required field, check if the field exists in the HDF5 file.
        for field in required_fields:
            check_h5py_File_has_group(h5file, groupname + field)
        
        ########################################################
        ### read data from HDF5 file and (optionally) close ####

        # output dictionary
        data = {
            "beta"        : None,
            "u_misc"      : None,
            "u_a"         : None,
            "u_d"         : None,
            "trait"       : None,
            "model_name"  : None,
            "hyperparams" : None,
        }

        ##################################
        ### read mandatory data fields ###

        # read beta array (ndarray dtype = any)
        data["beta"] = h5py_File_read_ndarray(h5file, groupname + "beta")

        # read u_misc array (ndarray dtype = any)
        data["u_misc"] = h5py_File_read_ndarray(h5file, groupname + "u_misc")

        # read u_a array (ndarray dtype = any)
        data["u_a"] = h5py_File_read_ndarray(h5file, groupname + "u_a")

        # read u_d array (ndarray dtype = any)
        data["u_d"] = h5py_File_read_ndarray(h5file, groupname + "u_d")

        #################################
        ### read optional data fields ###

        # read trat array (ndarray dtype = unicode / object)
        if groupname + "trait" in h5file:
            data["trait"] = h5py_File_read_ndarray_utf8(h5file, groupname + "trait")

        # read model_name data (dtype = str)
        if groupname + "model_name" in h5file:
            data["model_name"] = h5py_File_read_utf8(h5file, groupname + "model_name")

        # read hyperparams data (dtype = dict)
        if groupname + "hyperparams" in h5file:
            data["hyperparams"] = h5py_File_read_dict(h5file, groupname + "hyperparams")

        ######################
        ### close the file ###

        # close the file, only if the provided fieldname was a string or Path an not an h5py.File.
        if isinstance(filename, (str,Path)):
            h5file.close()

        ########################################################
        ################### Object creation ####################

        # create object from read data
        out = cls(
            beta        = data["beta"],
            u_misc      = data["u_misc"],
            u_a         = data["u_a"],
            u_d         = data["u_d"],
            trait       = data["trait"],
            model_name  = data["model_name"],
            hyperparams = data["hyperparams"],
        )

        return out



################################## Utilities ###################################
def check_is_DenseAdditiveDominanceLinearGenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type DenseAdditiveDominanceLinearGenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, DenseAdditiveDominanceLinearGenomicModel):
        raise TypeError("variable '{0}' must be a DenseAdditiveDominanceLinearGenomicModel".format(vname))
