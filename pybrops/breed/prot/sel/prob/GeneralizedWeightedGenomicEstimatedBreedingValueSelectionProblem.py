"""
Module implementing generalized weighted genomic selection as a subset optimization problem.
"""

__all__ = [
    "GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem",
    "GeneralizedWeightedGenomicEstimatedBreedingValueRealSelectionProblem",
    "GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelectionProblem"
]

from numbers import Integral, Number, Real
from typing import Callable, Optional, Union
import numpy
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SelectionProblem import SelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Real
from pybrops.core.error.error_value_numpy import check_ndarray_ndim
from pybrops.core.error.error_value_python import check_Number_in_interval

class GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem(SelectionProblem):
    """
    Helper class containing common properties for GWGS Problems.
    """
    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in wGEBV matrix
        return self._gwgebv.shape[1]

    @property
    def Z_a(self) -> numpy.ndarray:
        """Genotype matrix of shape (n,p)."""
        return self._Z_a
    @Z_a.setter
    def Z_a(self, value: numpy.ndarray) -> None:
        """Set genotype matrix."""
        check_is_ndarray(value, "Z")
        check_ndarray_ndim(value, "Z", 2)
        self._Z_a = value
    
    @property
    def u_a(self) -> numpy.ndarray:
        """Additive marker effects matrix of shape (p,t)."""
        return self._u_a
    @u_a.setter
    def u_a(self, value: numpy.ndarray) -> None:
        """Set additive marker effects matrix."""
        check_is_ndarray(value, "u_a")
        check_ndarray_ndim(value, "u_a", 2)
        self._u_a = value
    
    @property
    def fafreq(self) -> numpy.ndarray:
        """Favorable allele frequency matrix of shape (p,t)."""
        return self._fafreq
    @fafreq.setter
    def fafreq(self, value: numpy.ndarray) -> None:
        """Set favorable allele frequency matrix."""
        check_is_ndarray(value, "fafreq")
        check_ndarray_ndim(value, "fafreq", 2)
        # where there is a favorable allele frequency of 0,
        # convert to 1 to avoid division by zero
        value[value <= 0] = 1
        self._fafreq = value

    @property
    def alpha(self) -> Real:
        """Exponent to which to raise the favorable allele frequency. Must be in the range [0,1]."""
        return self._alpha
    @alpha.setter
    def alpha(self, value: Real) -> None:
        """Set exponent to which to raise the favorable allele frequency."""
        check_is_Real(value, "alpha")
        check_Number_in_interval(value, "alpha", 0, 1)
        self._alpha = value

    @property
    def gwgebv(self) -> numpy.ndarray:
        """Generalized weighted genomic estimated breeding values matrix of shape (n,t)."""
        return self._gwgebv
    @gwgebv.setter
    def gwgebv(self, value: numpy.ndarray) -> None:
        """Set generalized weighted genomic estimated breeding values matrix."""
        check_is_ndarray(value, "wgebv")
        check_ndarray_ndim(value, "wgebv", 2)
        self._gwgebv = value

class GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem(SubsetSelectionProblem,GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem):
    """
    docstring for SubsetGeneralizedWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            alpha: Real,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetGeneralizedWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(GeneralizedWeightedGenomicEstimatedBreedingValueSubsetSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )

        # set matrix values
        self.Z_a = Z_a
        self.u_a = u_a
        self.fafreq = fafreq
        self.alpha = alpha

        # calculate wGEBVs
        # (n,p) @ (p,t) -> (n,t)
        self.gwgebv = self.Z_a.dot(self.u_a * numpy.power(self.fafreq, -self.alpha))

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Generalized Weighted Genomic 
        Selection (GWGS). Scoring for GWGS is defined as the negative mean of 
        Weighted Genomic Estimated Breeding Values (wGEBV) for a selected subset.
        
        All objectives in this function are minimizing objectives; lower is better.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,)``.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A wGEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,t)[(ndecn,),:] -> (ndecn,t)    # select individuals
        # Step 2: (ndecn,t).sum(0)  -> (t,)         # sum across all individuals
        # Step 3: scalar * (t,) -> (t,)             # take mean across selection
        out = -(1.0 / len(x)) * (self._gwgebv[x,:].sum(0))

        return out

class GeneralizedWeightedGenomicEstimatedBreedingValueRealSelectionProblem(RealSelectionProblem,GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem):
    """
    docstring for RealGeneralizedWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            alpha: Real,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for RealGeneralizedWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(GeneralizedWeightedGenomicEstimatedBreedingValueRealSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )

        # set matrix values
        self.Z_a = Z_a
        self.u_a = u_a
        self.fafreq = fafreq
        self.alpha = alpha

        # calculate wGEBVs
        # (n,p) @ (p,t) -> (n,t)
        self.gwgebv = self.Z_a.dot(self.u_a * numpy.power(self.fafreq, -self.alpha))

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Generalized Weighted Genomic 
        Selection (GWGS). Scoring for GWGS is defined as the negative mean of 
        Weighted Genomic Estimated Breeding Values (wGEBV) for a population.

        All objectives in this function are minimizing objectives; lower is better.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A wGEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / x.sum()) * x

        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._gwgebv)

        return out

class GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelectionProblem(IntegerSelectionProblem,GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem):
    """
    docstring for IntegerGeneralizedWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            alpha: Real,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for IntegerGeneralizedWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(GeneralizedWeightedGenomicEstimatedBreedingValueIntegerSelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )

        # set matrix values
        self.Z_a = Z_a
        self.u_a = u_a
        self.fafreq = fafreq
        self.alpha = alpha

        # calculate wGEBVs
        # (n,p) @ (p,t) -> (n,t)
        self.gwgebv = self.Z_a.dot(self.u_a * numpy.power(self.fafreq, -self.alpha))

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Generalized Weighted Genomic 
        Selection (GWGS). Scoring for GWGS is defined as the negative mean of 
        Weighted Genomic Estimated Breeding Values (wGEBV) for a population.

        All objectives in this function are minimizing objectives; lower is better.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A wGEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / x.sum()) * x

        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._gwgebv)

        return out

class GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem(BinarySelectionProblem,GeneralizedWeightedGenomicEstimatedBreedingValueSelectionProblem):
    """
    docstring for BinaryGeneralizedWeightedGenomicSelectionProblem.
    """
    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            Z_a: numpy.ndarray,
            u_a: numpy.ndarray,
            fafreq: numpy.ndarray,
            alpha: Real,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Number,None],
            decn_space_upper: Union[numpy.ndarray,Number,None],
            nobj: Integral,
            obj_wt: Optional[Union[numpy.ndarray,Real]] = None,
            obj_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            obj_trans_kwargs: Optional[dict] = None,
            nineqcv: Optional[Integral] = None,
            ineqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            ineqcv_trans_kwargs: Optional[dict] = None,
            neqcv: Optional[Integral] = None,
            eqcv_wt: Optional[Union[numpy.ndarray,Real]] = None,
            eqcv_trans: Optional[Callable[[numpy.ndarray,numpy.ndarray,dict],numpy.ndarray]] = None,
            eqcv_trans_kwargs: Optional[dict] = None,
            **kwargs: dict
        ) -> None:
        """
        Constructor for BinaryGeneralizedWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(GeneralizedWeightedGenomicEstimatedBreedingValueBinarySelectionProblem, self).__init__(
            ndecn = ndecn,
            decn_space = decn_space,
            decn_space_lower = decn_space_lower,
            decn_space_upper = decn_space_upper,
            nobj = nobj,
            obj_wt = obj_wt,
            obj_trans = obj_trans,
            obj_trans_kwargs = obj_trans_kwargs,
            nineqcv = nineqcv,
            ineqcv_wt = ineqcv_wt,
            ineqcv_trans = ineqcv_trans,
            ineqcv_trans_kwargs = ineqcv_trans_kwargs,
            neqcv = neqcv,
            eqcv_wt = eqcv_wt,
            eqcv_trans = eqcv_trans,
            eqcv_trans_kwargs = eqcv_trans_kwargs,
            **kwargs
        )

        # set matrix values
        self.Z_a = Z_a
        self.u_a = u_a
        self.fafreq = fafreq
        self.alpha = alpha

        # calculate wGEBVs
        # (n,p) @ (p,t) -> (n,t)
        self.gwgebv = self.Z_a.dot(self.u_a * numpy.power(self.fafreq, -self.alpha))

    ############################################################################
    ############################## Object Methods ##############################
    ############################################################################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Generalized Weighted Genomic 
        Selection (GWGS). Scoring for GWGS is defined as the negative mean of 
        Weighted Genomic Estimated Breeding Values (wGEBV) for a population.

        All objectives in this function are minimizing objectives; lower is better.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
        args : tuple
            Additional non-keyword arguments.
        kwargs : dict
            Additional keyword arguments.
        
        Returns
        -------
        out : numpy.ndarray
            A wGEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / x.sum()) * x

        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._gwgebv)

        return out
