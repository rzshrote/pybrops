"""
Module defining optimization problems for binary optimal constribution selection.
"""

from numbers import Integral, Real
from typing import Callable, Optional, Union

import numpy
from pybrops.breed.prot.sel.prob.DenseSubsetSelectionProblem import DenseSubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_is_2d, check_ndarray_is_square
from pymoo.core.problem import ElementwiseEvaluationFunction, LoopedElementwiseEvaluation



class OptimalContributionSelectionProblem(DenseSubsetSelectionProblem):
    """
    Semi-abstract class representing an Optimal Contribution Selection Problem.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            bv: numpy.ndarray,
            C: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        Constructor for OptimalContributionSelectionProblem.
        
        Parameters
        ----------
        bv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        C : numpy.ndarray
            An upper triangle matrix of shape ``(n,n)`` resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
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
        # call DenseSubsetSelectionProblem constructor
        super(OptimalContributionSelectionProblem, self).__init__(
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
        # order dependent assignments
        self.bv = bv
        self.C = C

    ############################################################################
    ############################ Object Properties #############################
    ############################################################################
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in BV matrix plus 1
        return 1 + self._bv.shape[1]

    @property
    def bv(self) -> numpy.ndarray:
        """Breeding value matrix."""
        return self._bv
    @bv.setter
    def bv(self, value: numpy.ndarray) -> None:
        """Set breeding value matrix."""
        check_is_ndarray(value, "bv")
        check_ndarray_is_2d(value, "bv")
        self._bv = value
    @bv.deleter
    def bv(self) -> None:
        """Delete breeding value matrix."""
        del self._bv

    @property
    def C(self) -> numpy.ndarray:
        """Cholesky decomposition of the kinship matrix."""
        return self._C
    @C.setter
    def C(self, value: numpy.ndarray) -> None:
        """Set Cholesky decomposition of the kinship matrix."""
        check_is_ndarray(value, "C")
        check_ndarray_is_2d(value, "C")
        check_ndarray_is_square(value, "C")
        self._C = value
    @C.deleter
    def C(self) -> None:
        """Delete Cholesky decomposition of the kinship matrix."""
        del self._C

class SubsetOptimalContributionSelectionProblem(OptimalContributionSelectionProblem):
    """
    Class representing an Optimal Contribution Selection Problem for subset
    search spaces.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            bv: numpy.ndarray,
            C: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        Constructor for SubsetOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        bv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        C : numpy.ndarray
            An upper triangle matrix of shape ``(n,n)`` resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
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
        # call DenseSubsetSelectionProblem constructor
        super(SubsetOptimalContributionSelectionProblem, self).__init__(
            bv = bv,
            C = C,
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
        Encode a candidate solution for the given Problem into an ``l`` 
        dimensional latent evaluation space.
        
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
            A matrix of shape (1+t,).

            The first index in the array is the mean genomic relationship 
            (a minimizing objective):

            .. math::
                MGR = || \\textbf{C} \\textbf{(sel)} ||_2

            The next `t` indices in the array are the sum of breeding values for 
            each of ``t`` traits for the selection (all maximizing objectives).

            Where:

            - ``t`` is the number of traits.
        """
        # calculate MEH
        # (n,n)[:,(k,)] -> (n,k)
        # scalar * (n,k).sum(1) -> (n,)
        Cx = (1.0 / len(x)) * self.C[:,x].sum(1)

        # calculate mean genomic relationship
        # norm2( (n,), keepdims=True ) -> (1,)
        mgr = numpy.linalg.norm(Cx, ord = 2, keepdims = True)

        # calculate breeding value of the selection
        # (n,t)[(k,),:] -> (k,t)
        # (k,t).sum(0) -> (t,)
        gain = self.bv[x,:].sum(0)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgr,gain])

        # return (1+t,)
        return out

class RealOptimalContributionSelectionProblem(OptimalContributionSelectionProblem):
    """
    Class representing an Optimal Contribution Selection Problem for subset
    search spaces.
    """

    ############################################################################
    ########################## Special Object Methods ##########################
    ############################################################################
    def __init__(
            self,
            bv: numpy.ndarray,
            C: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        Constructor for RealOptimalContributionSelectionProblem.
        
        Parameters
        ----------
        bv : numpy.ndarray
            A breeding value matrix of shape ``(n,t)``. 
            If you are using a penalization transformation function, preferably
            these breeding values are centered and scaled to make the penalies 
            less extreme.

            Where:

            - ``n`` is the number of individuals.
            - ``t`` is the number of traits.
        C : numpy.ndarray
            An upper triangle matrix of shape ``(n,n)`` resulting from a Cholesky 
            decomposition of a kinship matrix: K = C'C.

            Where:

            - ``n`` is the number of individuals.
        ndecn : Integral
            Number of decision variables.
        decn_space: numpy.ndarray, None
            An array of shape ``(2,ndecn)`` defining the decision space.
            If None, do not set a decision space.
        decn_space_lower: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing lower limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
            If None, do not set a lower limit for the decision variables.
        decn_space_upper: numpy.ndarray, Real, None
            An array of shape ``(ndecn,)`` containing upper limits for decision variables.
            If a Real is provided, construct an array of shape ``(ndecn,)`` containing the Real.
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
        # call DenseSubsetSelectionProblem constructor
        super(RealOptimalContributionSelectionProblem, self).__init__(
            bv = bv,
            C = C,
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
        Encode a candidate solution for the given Problem into an ``l`` 
        dimensional latent evaluation space.
        
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
            A matrix of shape (1+t,).

            The first index in the array is the mean genomic relationship 
            (a minimizing objective):

            .. math::
                MGR = || \\textbf{C} \\textbf{(sel)} ||_2

            The next `t` indices in the array are the sum of breeding values for 
            each of ``t`` traits for the selection (all maximizing objectives).

            Where:

            - ``t`` is the number of traits.
        """
        # scale x to have a sum of 1 (contribution)
        contrib = (1.0 / x.sum()) * x

        # calculate mean genomic contribution
        # (n,n) . (n,) -> (n,)
        # scalar * (n,) -> (n,)
        # norm2( (n,), keepdims=True ) -> (1,)
        mgc = numpy.linalg.norm(self.C.dot(contrib), ord = 2, keepdims = True)

        # calculate breeding value of the selection
        # (n,t)[(k,),:] -> (k,t)
        # (k,t).sum(0) -> (t,)
        gain = self.bv.T.dot(contrib)
        
        # concatenate everything
        # (1,) concat (t,) -> (1+t,)
        out = numpy.concatenate([mgc,gain])

        # return (1+t,)
        return out

class IntegerOptimalContributionSelectionProblem(RealOptimalContributionSelectionProblem):
    # Real version is identical in implementation to Integer version
    pass



################################################################################
################################## Utilities ###################################
################################################################################
def check_is_OptimalContributionSelectionProblem(v: object, vname: str) -> None:
    """
    Check if object is of type OptimalContributionSelectionProblem, otherwise raise TypeError.

    Parameters
    ----------
    v : object
        Any Python object to test.
    vname : str
        Name of variable to print in TypeError message.
    """
    if not isinstance(v, OptimalContributionSelectionProblem):
        raise TypeError("'{0}' must be of type OptimalContributionSelectionProblem.".format(vname))