"""
Module defining basal interfaces and error checking routines for genomic models.
"""

from abc import ABCMeta, abstractmethod
from typing import Optional, Union
import numpy
from pybrops.core.io.HDF5InputOutput import HDF5InputOutput
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix
from pybrops.popgen.gmat.GenotypeMatrix import GenotypeMatrix
from pybrops.popgen.ptdf.PhenotypeDataFrame import PhenotypeDataFrame

class GenomicModel(HDF5InputOutput,metaclass=ABCMeta):
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

    ########################## Special Object Methods ##########################
    @abstractmethod
    def __copy__(
            self
        ) -> 'GenomicModel':
        """
        Make a shallow copy of the GenomicModel.

        Returns
        -------
        out : GenomicModel
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def __deepcopy__(
            self, 
            memo: dict
        ) -> 'GenomicModel':
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

    ############################ Object Properties #############################

    ################## Genomic Model Data ##################
    @property
    @abstractmethod
    def model_name(self) -> str:
        """Name of the model."""
        raise NotImplementedError("property is abstract")
    @model_name.setter
    @abstractmethod
    def model_name(self, value: str) -> None:
        """Set the name of the model"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def params(self) -> dict:
        """Model parameters."""
        raise NotImplementedError("property is abstract")
    @params.setter
    @abstractmethod
    def params(self, value: dict) -> None:
        """Set the model parameters"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def trait(self) -> numpy.ndarray:
        """Names of the traits predicted by the model."""
        raise NotImplementedError("property is abstract")
    @trait.setter
    @abstractmethod
    def trait(self, value: numpy.ndarray) -> None:
        """Set the names of the traits predicted by the model"""
        raise NotImplementedError("property is abstract")

    @property
    @abstractmethod
    def ntrait(self) -> int:
        """Number of traits predicted by the model."""
        raise NotImplementedError("property is abstract")
    @ntrait.setter
    @abstractmethod
    def ntrait(self, value: int) -> None:
        """Set the number of traits predicted by the model"""
        raise NotImplementedError("property is abstract")

    ############################## Object Methods ##############################

    ####### methods for model fitting and prediction #######
    @abstractmethod
    def fit_numpy(
            self, 
            Y: numpy.ndarray, 
            X: numpy.ndarray, 
            Z: numpy.ndarray, 
            **kwargs: dict
        ) -> None:
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

    @abstractmethod
    def fit(
            self, 
            ptobj: Union[BreedingValueMatrix,PhenotypeDataFrame,numpy.ndarray], 
            cvobj: numpy.ndarray, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            **kwargs: dict
        ) -> None:
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

    @abstractmethod
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
            A matrix of predicted breeding values.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def predict(
            self, 
            cvobj: numpy.ndarray, 
            gtobj: GenotypeMatrix, 
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

    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def score(
            self, 
            ptobj: Union[BreedingValueMatrix,PhenotypeDataFrame], 
            cvobj: object, 
            gtobj: GenotypeMatrix, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
    @abstractmethod
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def gebv(
            self, 
            gtobj: GenotypeMatrix, 
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
        gebvmat_hat : BreedingValueMatrix
            Genomic estimated breeding values.
        """
        raise NotImplementedError("method is abstract")

    ###### methods for population variance prediction ######
    @abstractmethod
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
            Array of shape ``(t,)`` contianing genetic variances for each trait.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def var_G(
            self, 
            gtobj: GenotypeMatrix, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
            Array of shape ``(t,)`` contianing genetic variances for each trait.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
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
            Array of shape ``(t,)`` contianing additive genetic variances for each trait.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def var_A(
            self, 
            gtobj: GenotypeMatrix, 
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
            Array of shape ``(t,)`` contianing additive genetic variances for each trait.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def var_a_numpy(
            self, 
            p: numpy.ndarray, 
            ploidy: int, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
            Array of shape ``(t,)`` contianing additive genic variances for each trait.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def var_a(
            self, 
            gtobj: Union[GenotypeMatrix,numpy.ndarray], 
            ploidy: int, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
            Array of shape ``(t,)`` contianing additive genic variances for each trait.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def bulmer_numpy(
            self, 
            Z: numpy.ndarray, 
            p: numpy.ndarray, 
            ploidy: int, 
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
            Array of shape ``(t,)`` contianing Bulmer effects for each trait.
            In the event that additive genic variance is zero, NaN's are produced.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def bulmer(
            self, 
            gtobj: GenotypeMatrix, 
            ploidy: int, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Calculate the Bulmer effect.

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
            Array of shape ``(t,)`` contianing Bulmer effects for each trait.
            In the event that additive genic variance is zero, NaN's are produced.
        """
        raise NotImplementedError("method is abstract")

    ############# methods for selection limits #############
    @abstractmethod
    def usl_numpy(
            self, 
            p: numpy.ndarray, 
            ploidy: int, 
            descale: bool, 
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
            An array of shape (t,) containing upper selection limits for each of ``t`` traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def usl(
            self, 
            gtobj: GenotypeMatrix, 
            ploidy: int, 
            descale: bool, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
            An array of shape (t,) containing upper selection limits for each of ``t`` traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def lsl_numpy(
            self, 
            p: numpy.ndarray, 
            ploidy: int, 
            descale: bool, 
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
            An array of shape (t,) containing lower selection limits for each of ``t`` traits.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def lsl(
            self, 
            gtobj: GenotypeMatrix, 
            ploidy: int, 
            descale: bool, 
            **kwargs: dict
        ) -> numpy.ndarray:
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
            An array of shape (t,) containing lower selection limits for each of ``t`` traits.
        """
        raise NotImplementedError("method is abstract")

    ############ methods for allele attributes #############
    @abstractmethod
    def facount(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype], 
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
            A numpy.ndarray of shape ``(p,)`` containing allele counts of the favorable allele.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def fafreq(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.ndarray], 
            **kwargs: dict
        ) -> numpy.ndarray:
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
            A numpy.ndarray of shape ``(p,t)`` containing allele frequencies of the favorable allele.
        """
        raise NotImplementedError("method is abstract")
    
    @abstractmethod
    def faavail(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.dtype] = None, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Determine whether a favorable allele is available in the present taxa.
        
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
            A numpy.ndarray of shape ``(p,t)`` containing whether a favorable 
            allele is available.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def fafixed(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.ndarray], 
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
            A numpy.ndarray of shape ``(p,)`` containing whether a favorable allele is fixed.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def dacount(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.ndarray], 
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
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def dafreq(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.ndarray], 
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
        raise NotImplementedError("method is abstract")
    
    @abstractmethod
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
            Genotype matrix for which to determine deleterious allele 
            frequencies.
        dtype : numpy.dtype, None
            Datatype of the returned array. If ``None``, use the native boolean 
            type.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A numpy.ndarray of shape ``(p,t)`` containing whether a deleterious 
            allele is available.
        """
        raise NotImplementedError("method is abstract")

    @abstractmethod
    def dafixed(
            self, 
            gmat: GenotypeMatrix, 
            dtype: Optional[numpy.ndarray], 
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
        raise NotImplementedError("method is abstract")



################################## Utilities ###################################
def check_is_GenomicModel(v: object, vname: str) -> None:
    """
    Check if object is of type GenomicModel. Otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, GenomicModel):
        raise TypeError("variable '{0}' must be a GenomicModel".format(vname))
