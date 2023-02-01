"""
Module defining basal interfaces and error checking routines for genomic models.
"""

from typing import Any, Union
import numpy
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix

class GenomicModel(HDF5InputOutput):
    """
    An abstract class for genomic models.

    The purpose for this abstract interface is to provide base functionality for:
        1) Model metadata storage.
        2) Model fitting.
        3) Model prediction.
        4) Model scoring.
        5) Genetic variance estimation using the model.
        6) Estimation of upper and lower selection limits using the model.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(self, **kwargs: dict) -> None:
        """
        Constructor for the GenomicModel class.

        Parameters
        ----------
        kwargs : dict
            Used for cooperative inheritance. Dictionary passing unused
            arguments to the parent class constructor.
        """
        super(GenomicModel, self).__init__(**kwargs)

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

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################

    ################## Genomic Model Data ##################
    def model_name():
        doc = "Name of the model"
        def fget(self):
            """Get the name of the model"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the name of the model"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the name of the model"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    model_name = property(**model_name())

    def params():
        doc = "Model parameters"
        def fget(self):
            """Get the model parameters"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the model parameters"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the model parameters"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    params = property(**params())

    def trait():
        doc = "Names of the traits predicted by the model"
        def fget(self):
            """Get the names of the traits predicted by the model"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the names of the traits predicted by the model"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the names of the traits predicted by the model"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    trait = property(**trait())

    def ntrait():
        doc = "Number of traits predicted by the model"
        def fget(self):
            """Get the number of traits predicted by the model"""
            raise NotImplementedError("method is abstract")
        def fset(self, value):
            """Set the number of traits predicted by the model"""
            raise NotImplementedError("method is abstract")
        def fdel(self):
            """Delete the number of traits predicted by the model"""
            raise NotImplementedError("method is abstract")
        return {"doc":doc, "fget":fget, "fset":fset, "fdel":fdel}
    ntrait = property(**ntrait())

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################

    ####### methods for model fitting and prediction #######
    def fit_numpy(self, Y, X, Z, **kwargs: dict):
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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def fit(self, ptobj, cvobj, gtobj, **kwargs: dict):
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
        kwargs : dict
            Additional keyword arguments.
        """
        raise NotImplementedError("method is abstract")

    def predict_numpy(self, X, Z, **kwargs: dict):
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
            A matrix of predicted breeding values.
        """
        raise NotImplementedError("method is abstract")

    def predict(self, cvobj, gtobj, **kwargs: dict):
        """
        Predict breeding values.

        Remark: The difference between 'predict' and 'gebv' is that 'predict'
        can incorporate other factors (e.g., fixed effects) to provide
        prediction estimates.

        Parameters
        ----------
        cvobj : numpy.ndarray
            An object containing covariate data.
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        bvmat_hat : BreedingValueMatrix
            Estimated breeding values.
        """
        raise NotImplementedError("method is abstract")

    def score_numpy(self, Y, X, Z, **kwargs: dict):
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
        raise NotImplementedError("method is abstract")

    def score(self, ptobj, cvobj, gtobj, **kwargs: dict):
        """
        Return the coefficient of determination R**2 of the prediction.

        Parameters
        ----------
        ptobj : BreedingValueMatrix or PhenotypeDataFrame
            An object containing phenotype data. Must be a matrix of breeding
            values or a phenotype data frame.
        cvobj : object
            An object containing covariate data.
        gtobj : GenotypeMatrix
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
        raise NotImplementedError("method is abstract")

    ######## methods for estimated breeding values #########
    def gebv_numpy(self, Z, **kwargs: dict):
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
        raise NotImplementedError("method is abstract")

    def gebv(self, gtobj, **kwargs: dict):
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
        gebvmat_hat : BreedingValueMatrix
            Genomic estimated breeding values.
        """
        raise NotImplementedError("method is abstract")

    ###### methods for population variance prediction ######
    def var_G_numpy(self, Z, **kwargs: dict):
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
        """
        raise NotImplementedError("method is abstract")

    def var_G(self, gtobj, **kwargs: dict):
        """
        Calculate the population genetic variance.

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
        """
        raise NotImplementedError("method is abstract")

    def var_A_numpy(self, Z, **kwargs: dict):
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
        """
        raise NotImplementedError("method is abstract")

    def var_A(self, gtobj, **kwargs: dict):
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
        """
        raise NotImplementedError("method is abstract")

    def var_a_numpy(self, p, ploidy, **kwargs: dict):
        """
        Calculate the population additive genic variance

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
        """
        raise NotImplementedError("method is abstract")

    def var_a(self, gtobj, ploidy, **kwargs: dict):
        """
        Calculate the population additive genic variance

        Parameters
        ----------
        gtobj : GenotypeMatrix, numpy.ndarray
            An object containing genotype data. Must be a matrix of genotype
            values.
        ploidy : int
            Ploidy of the species.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")

    def bulmer_numpy(self, Z, p, ploidy, **kwargs: dict):
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
        """
        raise NotImplementedError("method is abstract")

    def bulmer(self, gtobj, ploidy, **kwargs: dict):
        """
        Calculate the Bulmer effect.

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
            Array of Bulmer effects for each trait. In the event that additive
            genic variance is zero, NaN's are produced.
        """
        raise NotImplementedError("method is abstract")

    ############ methods for allele attributes #############
    def facount(self, gmat: GenotypeMatrix, dtype: Union[numpy.dtype,None], **kwargs: dict) -> numpy.ndarray:
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
            A numpy.ndarray of shape ``(p,)`` containing allele counts of the favorable allele.
        """
        raise NotImplementedError("method is abstract")

    def fafreq(self, gmat: GenotypeMatrix, dtype: Union[numpy.dtype,None], **kwargs: dict) -> numpy.ndarray:
        """
        Favorable allele frequency across all taxa.
        
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
            A numpy.ndarray of shape ``(p,)`` containing allele frequencies of the favorable allele.
        """
        raise NotImplementedError("method is abstract")
    
    def fafixed(self, gmat: GenotypeMatrix, dtype: Union[numpy.dtype,None], **kwargs: dict) -> numpy.ndarray:
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
            A numpy.ndarray of shape ``(p,)`` containing whether a favorable allele is fixed.
        """
        raise NotImplementedError("method is abstract")

    def dacount(self, gmat: GenotypeMatrix, dtype: Union[numpy.dtype,None], **kwargs: dict) -> numpy.ndarray:
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
        raise NotImplementedError("method is abstract")

    def dafreq(self, gmat: GenotypeMatrix, dtype: Union[numpy.dtype,None], **kwargs: dict) -> numpy.ndarray:
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
        raise NotImplementedError("method is abstract")
    
    def dafixed(self, gmat: GenotypeMatrix, dtype: Union[numpy.dtype,None], **kwargs: dict) -> numpy.ndarray:
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
        raise NotImplementedError("method is abstract")

    ############# methods for selection limits #############
    def usl_numpy(self, p, ploidy, descale, **kwargs: dict):
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
        """
        raise NotImplementedError("method is abstract")

    def usl(self, gtobj, ploidy, descale, **kwargs: dict):
        """
        Calculate the upper selection limit for a population.

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        ploidy : int
            Ploidy of the species.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")

    def lsl_numpy(self, p, ploidy, descale, **kwargs: dict):
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
        """
        raise NotImplementedError("method is abstract")

    def lsl(self, gtobj, ploidy, descale, **kwargs: dict):
        """
        Calculate the lower selection limit for a population.

        Parameters
        ----------
        gtobj : GenotypeMatrix
            An object containing genotype data. Must be a matrix of genotype
            values.
        ploidy : int
            Ploidy of the species.
        kwargs : dict
            Additional keyword arguments.

        Returns
        -------
        out : numpy.ndarray
        """
        raise NotImplementedError("method is abstract")



################################################################################
################################## Utilities ###################################
################################################################################
def is_GenomicModel(v: Any) -> bool:
    """
    Determine whether an object is a GenomicModel.

    Parameters
    ----------
    v : Any
        Any Python object to test.

    Returns
    -------
    out : bool
        True or False for whether v is a GenomicModel object instance.
    """
    return isinstance(v, GenomicModel)

def check_is_GenomicModel(v: Any, vname: str) -> None:
    """
    Check if object is of type GenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : Any
        Any Python object to test.
    varname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GenomicModel):
        raise TypeError("variable '{0}' must be a GenomicModel".format(vname))
