"""
Module implementing generalized weighted genomic selection as a subset optimization problem.
"""

from numbers import Integral, Number, Real
from typing import Callable, Optional, Union
import numpy
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_type_python import check_is_Real
from pybrops.core.error.error_value_numpy import check_ndarray_is_2d
from pybrops.core.error.error_value_python import check_Number_in_interval


class SubsetGeneralizedWeightedGenomicSelectionProblem(SubsetSelectionProblem):
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
            obj_wt: numpy.ndarray,
            obj_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            obj_trans_kwargs: Optional[dict],
            nineqcv: Integral,
            ineqcv_wt: numpy.ndarray,
            ineqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            ineqcv_trans_kwargs: Optional[dict],
            neqcv: Integral,
            eqcv_wt: numpy.ndarray,
            eqcv_trans: Optional[Callable[[numpy.ndarray,dict],numpy.ndarray]],
            eqcv_trans_kwargs: Optional[dict],
            **kwargs: dict
        ) -> None:
        """
        Constructor for SubsetGeneralizedWeightedGenomicSelectionProblem.
        
        Parameters
        ----------
        kwargs : dict
            Additional keyword arguments used for cooperative inheritance.
        """
        super(SubsetGeneralizedWeightedGenomicSelectionProblem, self).__init__(
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
            elementwise = True,
            elementwise_func = ElementwiseEvaluationFunction,
            elementwise_runner = LoopedElementwiseEvaluation(),
            **kwargs
        )

        # set matrix values
        self.Z_a = Z_a
        self.u_a = u_a
        self.fafreq = fafreq
        self.alpha = alpha

        # calculate wGEBVs
        # (n,p) @ (p,t) -> (n,t)
        self.wgebv = self.Z_a.dot(self.u_a * numpy.power(self.fafreq, -self.alpha))

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in wGEBV matrix
        return self._wgebv.shape[1]

    @property
    def Z_a(self) -> numpy.ndarray:
        """Genotype matrix of shape (n,p)."""
        return self._Z_a
    @Z_a.setter
    def Z_a(self, value: numpy.ndarray) -> None:
        """Set genotype matrix."""
        check_is_ndarray(value, "Z")
        check_ndarray_is_2d(value, "Z")
        self._Z_a = value
    @Z_a.deleter
    def Z_a(self) -> None:
        """Delete genotype matrix."""
        del self._Z_a
    
    @property
    def u_a(self) -> numpy.ndarray:
        """Additive marker effects matrix of shape (p,t)."""
        return self._u_a
    @u_a.setter
    def u_a(self, value: numpy.ndarray) -> None:
        """Set additive marker effects matrix."""
        check_is_ndarray(value, "u_a")
        check_ndarray_is_2d(value, "u_a")
        self._u_a = value
    @u_a.deleter
    def u_a(self) -> None:
        """Delete additive marker effects matrix."""
        del self._u_a
    
    @property
    def fafreq(self) -> numpy.ndarray:
        """Favorable allele frequency matrix of shape (p,t)."""
        return self._fafreq
    @fafreq.setter
    def fafreq(self, value: numpy.ndarray) -> None:
        """Set favorable allele frequency matrix."""
        check_is_ndarray(value, "fafreq")
        check_ndarray_is_2d(value, "fafreq")
        # where there is a favorable allele frequency of 0,
        # convert to 1 to avoid division by zero
        value[value <= 0] = 1
        self._fafreq = value
    @fafreq.deleter
    def fafreq(self) -> None:
        """Delete favorable allele frequency matrix."""
        del self._fafreq

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
    @alpha.deleter
    def alpha(self) -> None:
        """Delete exponent to which to raise the favorable allele frequency."""
        del self._alpha

    @property
    def wgebv(self) -> numpy.ndarray:
        """Weighted genomic estimated breeding values matrix of shape (n,t)."""
        return self._wgebv
    @wgebv.setter
    def wgebv(self, value: numpy.ndarray) -> None:
        """Set weighted genomic estimated breeding values matrix."""
        check_is_ndarray(value, "wgebv")
        check_ndarray_is_2d(value, "wgebv")
        self._wgebv = value
    @wgebv.deleter
    def wgebv(self) -> None:
        """Delete weighted genomic estimated breeding values matrix."""
        del self._wgebv

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
        Selection (GWGS). Scoring for GWGS is defined as the sum of Weighted
        Genomic Estimated Breeding Values (wGEBV) for a population.

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
        return self._wgebv[x,:].sum(0)
