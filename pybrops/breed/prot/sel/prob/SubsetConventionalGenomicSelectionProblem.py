"""
Module implementing conventional genomic selection as a subset optimization problem.
"""

from numbers import Integral, Number
from typing import Callable, Optional, Union

import numpy
from pybrops.breed.prot.sel.prob.DenseSubsetSelectionProblem import DenseSubsetSelectionProblem
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation

from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_2d


class SubsetConventionalGenomicSelectionProblem(DenseSubsetSelectionProblem):
    """
    docstring for SubsetConventionalGenomicSelectionProblem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            gebv: numpy.ndarray,
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
        Constructor for SubsetConventionalGenomicSelectionProblem.
        
        Parameters
        ----------
        gebv : numpy.ndarray
            An array of shape (n,t) containing genomic estimated breeding values.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Number, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Number is provided, construct an array of shape ``(ndecn,)`` containing the Number.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Number, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Number is provided, construct an array of shape ``(ndecn,)`` containing the Number.
            If None, do not set a upper limit for the decision variables.
        nobj: Integral
            Number of objectives.
        obj_wt: numpy.ndarray
            Objective function weights.
        obj_trans: Callable, None
            A transformation function transforming a latent space vector to an objective space vector.
            The transformation function must be of the form: ``obj_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the identity transformation function: copy the latent space vector to the objective space vector.
        obj_trans_kwargs: dict, None
            Keyword arguments for the latent space to objective space transformation function.
            If None, an empty dictionary is used.
        nineqcv: Integral,
            Number of inequality constraints.
        ineqcv_wt: numpy.ndarray,
            Inequality constraint violation weights.
        ineqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an inequality constraint violation vector.
            The transformation function must be of the form: ``ineqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        ineqcv_trans_kwargs: Optional[dict],
            Keyword arguments for the latent space to inequality constraint violation space transformation function.
            If None, an empty dictionary is used.
        neqcv: Integral
            Number of equality constraints.
        eqcv_wt: numpy.ndarray
            Equality constraint violation weights.
        eqcv_trans: Callable, None
            A transformation function transforming a latent space vector to an equality constraint violation vector.
            The transformation function must be of the form: ``eqcv_trans(x: numpy.ndarray, **kwargs) -> numpy.ndarray``
            If None, use the empty set transformation function: return an empty vector of length zero.
        eqcv_trans_kwargs: dict, None
            Keyword arguments for the latent space to equality constraint violation space transformation function.
            If None, an empty dictionary is used.
        kwargs : dict
            Additional keyword arguments passed to the parent class (DenseSubsetSelectionProblem) constructor.
        """
        super(SubsetConventionalGenomicSelectionProblem, self).__init__(
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
        # assignments
        self.gebv = gebv

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in GEBV matrix
        return self._gebv.shape[1]

    @property
    def gebv(self) -> numpy.ndarray:
        """Genomic estimated breeding values."""
        return self._gebv
    @gebv.setter
    def gebv(self, value: numpy.ndarray) -> None:
        """Set genomic estimated breeding values."""
        check_is_ndarray(value, "gebv")
        check_ndarray_is_2d(value, "gebv")
        self._gebv = value
    @gebv.deleter
    def gebv(self) -> None:
        """Delete genomic estimated breeding values."""
        del self._gebv

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
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the sum of
        Genomic Estimated Breeding Values (GEBV) for a population.

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
            A GEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,t)[(ndecn,),:] -> (ndecn,t)    # select individuals
        # Step 2: (ndecn,t).sum(0)  -> (t,)         # sum across all individuals
        return self._gebv[x,:].sum(0)
