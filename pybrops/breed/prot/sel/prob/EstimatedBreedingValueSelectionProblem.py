"""
Module implementing selection based on Estimated Breeding Values for multiple search space types.
"""

__all__ = [
    "EstimatedBreedingValueBinarySelectionProblem",
    "EstimatedBreedingValueIntegerSelectionProblem",
    "EstimatedBreedingValueRealSelectionProblem",
    "EstimatedBreedingValueSubsetSelectionProblem",
]

from abc import ABCMeta, abstractmethod
from numbers import Integral, Number, Real
from typing import Callable, Optional, Union

import numpy
from pybrops.breed.prot.sel.prob.BinarySelectionProblem import BinarySelectionProblem
from pybrops.breed.prot.sel.prob.IntegerSelectionProblem import IntegerSelectionProblem
from pybrops.breed.prot.sel.prob.RealSelectionProblem import RealSelectionProblem
from pybrops.breed.prot.sel.prob.SubsetSelectionProblem import SubsetSelectionProblem
from pybrops.core.error.error_type_numpy import check_is_ndarray
from pybrops.core.error.error_value_numpy import check_ndarray_axis_len_eq, check_ndarray_axis_len_gteq, check_ndarray_ndim
from pybrops.popgen.bvmat.BreedingValueMatrix import BreedingValueMatrix


class EstimatedBreedingValueSelectionProblemMixin(
        metaclass = ABCMeta,
    ):
    """Helper class containing common properties for EBV selection problems."""

    ########################## Special Object Methods ##########################
    # __init__() CANNOT be defined to be classified as a Mixin class

    ############################ Object Properties #############################
    
    ############## Number of latent variables ##############
    @property
    def nlatent(self) -> Integral:
        """Number of latent variables."""
        # return number of traits in GEBV matrix
        return self._ebv.shape[1]

    ##################### GEBV matrix ######################
    @property
    def ebv(self) -> numpy.ndarray:
        """Genomic estimated breeding values."""
        return self._ebv
    @ebv.setter
    def ebv(self, value: numpy.ndarray) -> None:
        """Set genomic estimated breeding values."""
        check_is_ndarray(value, "ebv")
        check_ndarray_ndim(value, "ebv", 2)
        # most (binary, real, integer) problems require decisons for each cross
        check_ndarray_axis_len_eq(value, "ebv", 0, self.ndecn)
        self._ebv = value

    ############################## Class Methods ###############################
    @classmethod
    @abstractmethod
    def from_bvmat(
            cls,
            bvmat: BreedingValueMatrix,
            unscale: bool,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        ) -> "EstimatedBreedingValueSelectionProblemMixin":
        raise NotImplementedError("class method is abstract")

class EstimatedBreedingValueSubsetSelectionProblem(EstimatedBreedingValueSelectionProblemMixin,SubsetSelectionProblem):
    """
    Class representing selection on Estimated Breeding Values (EBVs) in subset search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        Constructor for EstimatedBreedingValueSubsetSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            An array of shape (n,t) containing estimated breeding values.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(EstimatedBreedingValueSubsetSelectionProblem, self).__init__(
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
        # assignments
        self.ebv = ebv

    ############################ Object Properties #############################
    
    ##################### GEBV matrix ######################
    @EstimatedBreedingValueSelectionProblemMixin.ebv.setter
    def ebv(self, value: numpy.ndarray) -> None:
        """Set genomic estimated breeding values."""
        check_is_ndarray(value, "ebv")
        check_ndarray_ndim(value, "ebv", 2)
        # for subset problems, must have more crosses than decision variables
        check_ndarray_axis_len_gteq(value, "ebv", 0, self.ndecn)
        self._ebv = value

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the mean of
        Genomic Estimated Breeding Values (GEBV) for a selected subset.

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
        # Step 3: scalar * (t,) -> (t,)             # take mean across selection
        out = -(1.0 / len(x)) * (self._ebv[x,:].sum(0))

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_bvmat(
            cls,
            bvmat: BreedingValueMatrix,
            unscale: bool,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        ) -> "EstimatedBreedingValueSubsetSelectionProblem":
        # extract EBVs for all individuals
        ebv = bvmat.unscale() if unscale else bvmat.mat

        # construct class
        out = cls(
            ebv = ebv,
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

        return out

class EstimatedBreedingValueRealSelectionProblem(EstimatedBreedingValueSelectionProblemMixin,RealSelectionProblem):
    """
    Class representing selection on Estimated Breeding Values (EBVs) in real search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for EstimatedBreedingValueRealSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            An array of shape (n,t) containing estimated breeding values.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(EstimatedBreedingValueRealSelectionProblem, self).__init__(
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
        # assignments
        self.ebv = ebv

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the mean of
        Genomic Estimated Breeding Values (GEBV) for a set of selection 
        contributions.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(k,) == (ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(ax)`` where ``a`` is any number.
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
        # calculate sum(x)
        xsum = x.sum()

        # if sum(x) ~== 0, then set to 1
        xsum = xsum if abs(xsum) >= 1e-10 else 1.0

        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / xsum) * x

        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._ebv)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_bvmat(
            cls,
            bvmat: BreedingValueMatrix,
            unscale: bool,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        ) -> "EstimatedBreedingValueRealSelectionProblem":
        # extract EBVs for all individuals
        ebv = bvmat.unscale() if unscale else bvmat.mat

        # construct class
        out = cls(
            ebv = ebv,
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

        return out

class EstimatedBreedingValueIntegerSelectionProblem(EstimatedBreedingValueSelectionProblemMixin,IntegerSelectionProblem):
    """
    Class representing selection on Estimated Breeding Values (EBVs) in integer search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for EstimatedBreedingValueIntegerSelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            An array of shape (n,t) containing estimated breeding values.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(EstimatedBreedingValueIntegerSelectionProblem, self).__init__(
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
        # assignments
        self.ebv = ebv

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the mean of
        Genomic Estimated Breeding Values (GEBV) for a set of selection 
        contributions.

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
            A GEBV matrix of shape ``(t,)``.

            Where:

            - ``t`` is the number of traits.
        """
        # calculate sum(x)
        xsum = x.sum()

        # if sum(x) ~== 0, then set to 1
        xsum = xsum if abs(xsum) >= 1e-10 else 1.0

        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / xsum) * x

        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._ebv)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_bvmat(
            cls,
            bvmat: BreedingValueMatrix,
            unscale: bool,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        ) -> "EstimatedBreedingValueIntegerSelectionProblem":
        # extract EBVs for all individuals
        ebv = bvmat.unscale() if unscale else bvmat.mat

        # construct class
        out = cls(
            ebv = ebv,
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

        return out

class EstimatedBreedingValueBinarySelectionProblem(EstimatedBreedingValueSelectionProblemMixin,BinarySelectionProblem):
    """
    Class representing selection on Estimated Breeding Values (EBVs) in binary search spaces.
    """

    ########################## Special Object Methods ##########################
    def __init__(
            self,
            ebv: numpy.ndarray,
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
        Constructor for EstimatedBreedingValueBinarySelectionProblem.
        
        Parameters
        ----------
        ebv : numpy.ndarray
            An array of shape (n,t) containing estimated breeding values.
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
            Additional keyword arguments passed to the parent class (SubsetSelectionProblem) constructor.
        """
        super(EstimatedBreedingValueBinarySelectionProblem, self).__init__(
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
        # assignments
        self.ebv = ebv

    ############################## Object Methods ##############################
    def latentfn(
            self, 
            x: numpy.ndarray, 
            *args: tuple, 
            **kwargs: dict
        ) -> numpy.ndarray:
        """
        Score a population of individuals based on Conventional Genomic Selection
        (CGS) (Meuwissen et al., 2001). Scoring for CGS is defined as the mean of
        Genomic Estimated Breeding Values (GEBV) for a set of selection 
        contributions.

        Parameters
        ----------
        x : numpy.ndarray
            A candidate solution vector of shape ``(ndecn,) == (ntaxa,)``.
            On entry, this vector is scaled to have a unit sum, such that
            ``latentfn(x) == latentfn(kx)`` where ``k`` is any number.
            In the case that ``sum(x)`` is approximately zero, then do not scale
            ``x`` to avoid division by zero.
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
        # calculate sum(x)
        xsum = x.sum()

        # if sum(x) ~== 0, then set to 1
        xsum = xsum if abs(xsum) >= 1e-10 else 1.0

        # scale x to have a sum of 1 (contribution)
        # (n,) -> (n,)
        contrib = (1.0 / xsum) * x

        # select individuals and take the sum of their GEBVs
        # CGS calculation explanation
        # Step 1: (n,) . (n,t) -> (t,)  # take dot product with contributions
        out = -contrib.dot(self._ebv)

        return out

    ############################## Class Methods ###############################
    @classmethod
    def from_bvmat(
            cls,
            bvmat: BreedingValueMatrix,
            unscale: bool,
            ndecn: Integral,
            decn_space: Union[numpy.ndarray,None],
            decn_space_lower: Union[numpy.ndarray,Real,None],
            decn_space_upper: Union[numpy.ndarray,Real,None],
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
        ) -> "EstimatedBreedingValueBinarySelectionProblem":
        # extract EBVs for all individuals
        ebv = bvmat.unscale() if unscale else bvmat.mat

        # construct class
        out = cls(
            ebv = ebv,
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

        return out
